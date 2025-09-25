import os
import subprocess
import pandas as pd
import numpy as np

def process_clustered_contigs(clstr_file, diamond_tsv, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    # Run clstr2txt.pl and capture to text file
    result = subprocess.run(
        ["clstr2txt.pl", clstr_file],
        capture_output=True,
        text=True,
        check=True
    )
    cluster_txt_path = os.path.join(output_dir, "clusters.txt")
    with open(cluster_txt_path, "w") as f:
        f.write(result.stdout)

    # --- Helpers -------------------------------------------------------------
    def pct_to_float(series):
        """Convert values like 97 or '97%' (or NaN) to float 97.0."""
        s = series.astype(str).str.replace('%', '', regex=False)
        # empty strings -> NaN
        s = s.replace({'': np.nan, 'nan': np.nan, 'None': np.nan})
        return pd.to_numeric(s, errors='coerce')

    def to_int01(series):
        """Convert representative flag to integer 0/1 robustly."""
        mapping = {'1': 1, '0': 0, 'True': 1, 'False': 0, 'Y': 1, 'N': 0, '*': 1}
        s = series.astype(str).str.strip().map(mapping).astype('float')
        return s.fillna(0).astype(int)

    def extract_len_from_query_id(qs: pd.Series):
        """
        Try to extract contig length from query_id.
        Works with patterns like:
          sample|contig_xxx_len_12345
        and falls back to the 4th underscore token (index 3) like your previous code.
        """
        qs = qs.astype(str)
        tail = qs.str.split('|').str[-1]
        # First try: regex 'len_12345' / 'len-12345'
        rx = tail.str.extract(r'len[_-]?(\d+)', expand=False)
        out = pd.to_numeric(rx, errors='coerce')
        # Fallback: 4th underscore-separated token
        fallback = pd.to_numeric(tail.str.split('_').str[3], errors='coerce')
        return out.fillna(fallback)

    # --- Load and normalize cluster table -----------------------------------
    # Read everything as string first to avoid .str accessor failures
    df = pd.read_csv(cluster_txt_path, sep="\t", dtype=str)

    # Safely parse percent columns (present in typical clstr2txt output)
    if "clstr_iden" in df.columns:
        df["clstr_iden"] = pct_to_float(df["clstr_iden"])
    if "clstr_cov" in df.columns:
        df["clstr_cov"] = pct_to_float(df["clstr_cov"])

    # Representative flag to 0/1
    if "clstr_rep" in df.columns:
        df["clstr_rep"] = to_int01(df["clstr_rep"])
    else:
        # If missing, default everyone to member (0)
        df["clstr_rep"] = 0

    # Select representative and high-coverage members
    # If clstr_cov missing, keep all non-rep members
    if "clstr_cov" in df.columns:
        high_cov_members = df[(df["clstr_rep"] == 0) & (df["clstr_cov"] >= 10)]
    else:
        high_cov_members = df[df["clstr_rep"] == 0]

    representatives = df[df["clstr_rep"] == 1]
    filtered = pd.concat([representatives, high_cov_members], ignore_index=True)
    if "clstr" in filtered.columns and "clstr_rep" in filtered.columns:
        filtered = filtered.sort_values(by=["clstr", "clstr_rep"], ascending=[True, False])

    # Extract metadata from original IDs
    if "id" not in filtered.columns:
        raise ValueError("Expected column 'id' from clstr2txt.pl output was not found.")
    filtered['Sample_ID'] = filtered['id'].astype(str).str.split('|').str[0]
    # Derive contig length from ID pattern
    filtered['contigs_len'] = extract_len_from_query_id(filtered['id'])

    # --- Load Diamond results and ensure qcov --------------------------------
    diamond = pd.read_csv(diamond_tsv, sep="\t", dtype=str)

    # Make sure we have contig length in Diamond df if we need to compute qcov
    if "contigs_len" not in diamond.columns and "query_id" in diamond.columns:
        diamond["contigs_len"] = extract_len_from_query_id(diamond["query_id"])

    # Numeric helpers for Diamond fields we may use
    for col in ("qstart", "qend", "aln_len", "contigs_len"):
        if col in diamond.columns:
            diamond[col] = pd.to_numeric(diamond[col], errors='coerce')

    # Compute qcov if missing:
    if "qcov" not in diamond.columns:
        qcov = None
        if {"qstart", "qend", "contigs_len"}.issubset(diamond.columns):
            qcov = ((diamond["qend"] - diamond["qstart"] + 1) / diamond["contigs_len"]) * 100
        elif {"aln_len", "contigs_len"}.issubset(diamond.columns):
            qcov = (diamond["aln_len"] / diamond["contigs_len"]) * 100
        elif {"aln_len", "qlen"}.issubset(diamond.columns):
            # if DIAMOND was run with qlen in outfmt
            diamond["qlen"] = pd.to_numeric(diamond["qlen"], errors='coerce')
            qcov = (diamond["aln_len"] / diamond["qlen"]) * 100

        if qcov is not None:
            diamond["qcov"] = qcov

    # Columns to pull from diamond results (take only what exists)
    diamond_cols = [
        'Sample_ID', 'query_id', 'pident', 'contigs_len', 'virus', 'evalue', 'bitscore',
        'aln_len', 'mismatches', 'gaps', 'qstart', 'qend', 'sstart', 'send', 'qcov'
    ]
    available_cols = [c for c in diamond_cols if c in diamond.columns]

    # Merge annotations from representative contigs onto their cluster members
    rep_annot = representatives.rename(columns={"id": "query_id"})
    rep_annot = pd.merge(rep_annot, diamond[available_cols], on='query_id', how='left')

    # Build mapping dicts keyed by cluster id
    if "clstr" not in rep_annot.columns:
        raise ValueError("Expected 'clstr' column not present in representative annotations.")
    cluster_maps = {
        col: dict(zip(rep_annot['clstr'], rep_annot[col]))
        for col in available_cols if col not in ['query_id', 'Sample_ID']
    }

    filtered = filtered.rename(columns={"id": "query_id"})
    for col, cmap in cluster_maps.items():
        filtered[col] = filtered['clstr'].map(cmap)

    # Keep only viral contigs if 'virus' column exists; otherwise keep all
    if 'virus' in filtered.columns:
        filtered = filtered[filtered['virus'].astype(str).str.contains("virus", case=False, na=False)]

    # Save
    output_file = os.path.join(output_dir, "filtered_clusters_assigned_rep_virus.tsv")
    filtered.to_csv(output_file, sep="\t", index=False)

    return output_file
