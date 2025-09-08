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

    df = pd.read_csv(cluster_txt_path, sep="\t")
    df["clstr_iden"] = df["clstr_iden"].str.rstrip('%').astype(float)
    df["clstr_cov"] = df["clstr_cov"].str.rstrip('%').astype(float)

    # Select representative and high-coverage members
    representatives = df[df["clstr_rep"] == 1]
    high_cov_members = df[(df["clstr_rep"] == 0) & (df["clstr_cov"] >= 70)]
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
    if 'qcov' not in diamond.columns:
        diamond['qcov'] = ((diamond['qend'] - diamond['qstart'] + 1) / diamond['contigs_len']) * 100

    # Merge representatives with Diamond to get annotations
    rep_annot = representatives.rename(columns={"id": "query_id"})
    rep_annot = pd.merge(
        rep_annot,
        diamond[['Sample_ID', 'query_id', 'pident', 'contigs_len', 'virus', 'evalue', 'bitscore',
                 'aln_len', 'mismatches', 'gaps', 'qstart', 'qend', 'sstart', 'send', 'qcov']],
        on='query_id', how='left'
    )

    # Create cluster-to-virus, evalue, bitscore, qcov mappings
    cluster_to_virus = dict(zip(rep_annot['clstr'], rep_annot['virus']))
    cluster_to_evalue = dict(zip(rep_annot['clstr'], rep_annot['evalue']))
    cluster_to_bitscore = dict(zip(rep_annot['clstr'], rep_annot['bitscore']))
    cluster_to_id = dict(zip(rep_annot['clstr'], rep_annot['pident']))
    cluster_to_qcov = dict(zip(rep_annot['clstr'], rep_annot['qcov']))

    # Map annotations back to filtered contigs
    filtered['virus'] = filtered['clstr'].map(cluster_to_virus)
    filtered['evalue'] = filtered['clstr'].map(cluster_to_evalue)
    filtered['bitscore'] = filtered['clstr'].map(cluster_to_bitscore)
    filtered['pident'] = filtered['clstr'].map(cluster_to_id)
    filtered['qcov'] = filtered['clstr'].map(cluster_to_qcov)

    # Filter to viral contigs with qcov ≥ 70%
    filtered = filtered.rename(columns={"id": "query_id"})
    filtered = filtered[filtered['virus'].str.contains("virus", case=False, na=False)]
    #filtered = filtered[
    #filtered['virus'].str.contains("virus", case=False, na=False) &
      #  (filtered['qcov'] >= 70)
   # ]

    # Output
    filtered_output_file = os.path.join(output_dir, "filtered_clusters_assigned_rep_virus.tsv")
    filtered.to_csv(filtered_output_file, sep="\t", index=False)

    return filtered_output_file
