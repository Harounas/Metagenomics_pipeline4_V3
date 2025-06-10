#!/usr/bin/env python3
"""
extract_contigs_diamond.py

This script extracts viral contigs from Kraken2 results, merges and renames them,
runs Diamond BLASTX against a specified database, and annotates hits via NCBI Entrez.
"""

from Bio import SeqIO, Entrez
import csv
import os
from pathlib import Path
import subprocess
import pandas as pd

# Set your email for NCBI Entrez usage
Entrez.email = "harounasoum17@gmail.com"

from pathlib import Path
from Bio import SeqIO
import csv

def extract_and_merge_contigs(base_contigs_dir, output_fasta="merged_contigs.fasta", min_length=200):
    """
    Extracts contigs >200 bp from each sample's contigs.fasta,
    regardless of taxon, and merges them into one output FASTA.

    Expects:
      - contigs.fasta files in {base_contigs_dir}/{sample_id}/
    """
    base_dir = Path(base_contigs_dir)
    merged_records = []

    for contigs_fasta in base_dir.glob("*/contigs.fasta"):
        sid = contigs_fasta.parent.name
        print(f"Processing {sid}")

        for rec in SeqIO.parse(str(contigs_fasta), "fasta"):
            if len(rec.seq) > min_length:
                rec.id = f"{sid}|{rec.id}"  # Prefix with sample ID for clarity
                rec.description = ""  # Remove long description
                merged_records.append(rec)

    if merged_records:
        SeqIO.write(merged_records, output_fasta, "fasta")
        print(f"✅ Merged {len(merged_records)} contigs >{min_length} bp into {output_fasta}")
    else:
        print("⚠️ No contigs found >200 bp in any sample.")

# Example usage
# extract_and_merge_contigs("/path/to/base_contigs_dir")

import subprocess
from pathlib import Path

def run_genomad_and_cluster(input_fasta, output_dir, genomad_db, final_output="clustered_contigs.fasta", 
                            min_length=200, identity=0.95, word_size=10, 
                            mem_mb=16000, threads=8):
    """
    Runs geNomad end-to-end and clusters viral contigs using cd-hit-est.
    
    Args:
        input_fasta (str): Path to input merged contigs FASTA file.
        output_dir (str): Path to store geNomad output.
        genomad_db (str): Path to geNomad database directory.
        final_output (str): Path for the final clustered FASTA output.
        min_length (int): Minimum contig length for geNomad.
        identity (float): Identity threshold for cd-hit-est.
        word_size (int): Word size for cd-hit-est.
        mem_mb (int): Memory (MB) for cd-hit-est.
        threads (int): Number of threads.
    """
    # 1️⃣ geNomad end-to-end
    genomad_cmd = [
        "genomad", "end-to-end",
        input_fasta, output_dir,
        "--db-dir", genomad_db,
        "--min-length", str(min_length),
        "--threads", str(threads)
    ]
    print(f"🚀 Running geNomad:\n{' '.join(genomad_cmd)}")
    subprocess.run(genomad_cmd, check=True)

    # 2️⃣ Path to geNomad's viral contigs output
    virus_fasta = Path(output_dir) /"merged_contigs_filtered_summary/merged_contigs_virus.fna"
    if not virus_fasta.exists():
        raise FileNotFoundError(f"Viral contigs FASTA not found: {virus_fasta}")

    # 3️⃣ cd-hit-est clustering
    cdhit_cmd = [
        "cd-hit-est",
        "-i", str(virus_fasta),
        "-o", output_dir/final_output,
        "-c", str(identity),
        "-n", str(word_size),
        "-d", "0",
        "-M", str(mem_mb),
        "-T", str(threads)
    ]
    print(f"🚀 Running cd-hit-est:\n{' '.join(cdhit_cmd)}")
    subprocess.run(cdhit_cmd, check=True)

    print(f"\n✅ Pipeline complete! Final clustered FASTA: {output_dir/final_output}")

# Example usage
# run_genomad_and_cluster(
#     "../../250328_VH01476_21_AACN7NGHV/merged_filtered_contigs.fasta",
#     "output_dir",
#     "../genomad_db/",
#     "output_dir/clustered_contigs.fna"
# )

"""
def merge_and_rename_contigs(base_contigs_dir, merged_filename=f"{output_dir}/{final_output}"):
    base_dir = Path(base_contigs_dir)
    merged = Path(merged_filename)

    with merged.open("w") as out:
        for sample_dir in base_dir.iterdir():
            if not sample_dir.is_dir():
                continue
            sid = sample_dir.name
            for fp in sample_dir.glob(f"{sid}_*.fasta"):
                with fp.open() as inf:
                    for line in inf:
                        if line.startswith(">"):
                            hdr = line[1:].strip()
                            out.write(f">{sid}|{hdr}\n")
                        else:
                            out.write(line)
                print(f"  Merged {fp}")
    print(f"\nMerged FASTA → {merged_filename}")
"""
def run_diamond(diamond_db, query_file="merged_contigs_renamed.fasta",
                output_file="results.m8", threads=8):
    cmd = [
        "diamond", "blastx",
        "--query", query_file,
        "--db", diamond_db,
        "--out", output_file,
        "--threads", str(threads)
    ]
    print("\nExecuting:", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True)
        print(f"Diamond results → {output_file}")
    except subprocess.CalledProcessError as e:
        print("Diamond failed:", e)

def process_diamond_results(results_filename="results.m8",
                            extracted_csv="extracted_virus.csv",
                            extracted_csv1="extracted_virus1.csv"):
    df = pd.read_csv(results_filename, sep="\t", header=None)
    df.columns = [
        'query_id','subject_id','pident','aln_len',
        'mismatches','gaps','qstart','qend','sstart','send',
        'evalue','bitscore'
    ]

    def fetch_name(acc):
        try:
            h = Entrez.efetch(db="protein", id=acc, rettype="gb", retmode="text")
            txt = h.read(); h.close()
            for l in txt.splitlines():
                if l.startswith("  ORGANISM"):
                    return l.split("ORGANISM")[-1].strip()
        except:
            return "Unknown"
        return "Unknown"

    df["Virus Name"] = df["subject_id"].apply(fetch_name)
    df.to_csv(extracted_csv, index=False)
    df['Custom_Score'] = (df.pident * df.aln_len) / (1 + df.mismatches + df.gaps)
    df = df.sort_values("Custom_Score", ascending=False)
    df.to_csv(extracted_csv1, index=False)
    print(f"\nAnnotated Diamond → {extracted_csv}, sorted → {extracted_csv1}")
