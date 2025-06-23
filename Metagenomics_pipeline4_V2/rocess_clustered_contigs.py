import pandas as pd
import subprocess
import os

def process_clustered_contigs(clstr_file, diamond_tsv, output_dir):
    """
    Processes .clstr file, filters contigs by coverage, annotates them with representative viruses,
    and outputs a filtered TSV file for downstream alignment.

    Parameters:
    - clstr_file (str): Path to the .clstr file (e.g., clustered_contigs.fasta.clstr)
    - diamond_tsv (str): Path to the Diamond result file with sample and virus info
    - output_dir (str): Directory where intermediate and final outputs will be written

    Returns:
    - filtered_output_file (str): Path to the final filtered TSV file
    """
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
    high_cov_members = df[(df["clstr_rep"] == 0) & (df["clstr_cov"] >= 10)]
    filtered = pd.concat([representatives, high_cov_members]).sort_values(
        by=["clstr", "clstr_rep"], ascending=[True, False]
    )

    # Extract metadata
    filtered['Sample_ID'] = filtered['id'].str.split('|').str[0]
    contig_part = filtered['id'].str.split('|').str[-1]
    filtered['contigs_len'] = pd.to_numeric(contig_part.str.split('_').str[3], errors='coerce')

    # Load Diamond results and merge virus annotations for representatives
    diamond = pd.read_csv(diamond_tsv, sep="\t")
    rep_virus = representatives.rename(columns={"id": "query_id"})
    rep_virus = pd.merge(rep_virus, diamond[['Sample_ID', 'query_id', 'contigs_len', 'virus']],
                         on='query_id', how='left')

    cluster_to_virus = dict(zip(rep_virus['clstr'], rep_virus['virus']))
    filtered['virus'] = filtered['clstr'].map(cluster_to_virus)

    # Filter for viral contigs and export
    filtered = filtered.rename(columns={"id": "query_id"})
    filtered = filtered[filtered['virus'].str.contains("virus", case=False, na=False)]

    filtered_output_file = os.path.join(output_dir, "filtered_clusters_assigned_rep_virus.tsv")
    filtered.to_csv(filtered_output_file, sep="\t", index=False)

    return filtered_output_file
