import pandas as pd
import subprocess
import os

def process_clustered_contigs(clstr_file, diamond_tsv, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    # Run clstr2txt.pl
    result = subprocess.run(
        ["clstr2txt.pl", clstr_file],
        capture_output=True,
        text=True,
        check=True
    )

    cluster_txt_path = os.path.join(output_dir, "clusters.txt")
    with open(cluster_txt_path, "w") as f:
        f.write(result.stdout)

    # Load cluster data
    df = pd.read_csv(cluster_txt_path, sep="\t")
    df["clstr_iden"] = df["clstr_iden"].str.rstrip('%').astype(float)
    df["clstr_cov"] = df["clstr_cov"].str.rstrip('%').astype(float)

    # Select representative and high-coverage members
    representatives = df[df["clstr_rep"] == 1]
    high_cov_members = df[(df["clstr_rep"] == 0) & (df["clstr_cov"] >= 10)]
    filtered = pd.concat([representatives, high_cov_members]).sort_values(
        by=["clstr", "clstr_rep"], ascending=[True, False]
    )

    # Extract metadata
    filtered['Sample_ID'] = filtered['id'].str.split('|').str[0]
    contig_part = filtered['id'].str.split('|').str[-1]
    filtered['contigs_len'] = pd.to_numeric(contig_part.str.split('_').str[3], errors='coerce')

    # Load Diamond results
    diamond = pd.read_csv(diamond_tsv, sep="\t")

    # Calculate query coverage if not provided
    if 'qcov' not in diamond.columns and 'qend' in diamond.columns and 'qstart' in diamond.columns:
        diamond['qcov'] = ((diamond['qend'] - diamond['qstart'] + 1) / diamond['contigs_len']) * 100

    # Columns to pull from diamond results
    diamond_cols = [
        'Sample_ID', 'query_id', 'pident', 'contigs_len', 'virus', 'evalue', 'bitscore',
        'aln_len', 'mismatches', 'gaps', 'qstart', 'qend', 'sstart', 'send', 'qcov'
    ]
    available_cols = [c for c in diamond_cols if c in diamond.columns]

    # Merge annotations to representative contigs
    rep_annot = representatives.rename(columns={"id": "query_id"})
    rep_annot = pd.merge(rep_annot, diamond[available_cols], on='query_id', how='left')

    # Build mapping dictionaries for each annotation field
    cluster_maps = {
        col: dict(zip(rep_annot['clstr'], rep_annot[col])) 
        for col in available_cols if col not in ['query_id', 'Sample_ID']
    }

    # Apply annotation maps to all filtered cluster members
    filtered = filtered.rename(columns={"id": "query_id"})
    for col, cmap in cluster_maps.items():
        filtered[col] = filtered['clstr'].map(cmap)

    # Keep only viral contigs
    filtered = filtered[filtered['virus'].str.contains("virus", case=False, na=False)]

    # Save output
    output_file = os.path.join(output_dir, "filtered_clusters_assigned_rep_virus.tsv")
    filtered.to_csv(output_file, sep="\t", index=False)

    return output_file
