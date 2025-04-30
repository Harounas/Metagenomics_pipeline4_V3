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

def extract_contigs(base_contigs_dir, summary_filename="contigs_summary.tsv"):
    """
    Extract viral contigs (>500 bp) based on Kraken2 reports, write per-taxon FASTA files,
    and save a summary TSV.

    Expects:
      - {sample_id}_Viruses_kraken_report.txt (in base_contigs_dir)
      - {sample_id}_kraken2_output.txt (in base_contigs_dir)
      - contigs.fasta (in base_contigs_dir/{sample_id}/)
    """
    allowed_ranks = {"F", "F1", "F2", "G", "G1", "G2", "S", "S1", "S2"}
    base_dir = Path(base_contigs_dir)

    with open(summary_filename, "w") as summary_f:
        summary_f.write("sample_id\ttaxon_id\tpathogen\tlong_contigs\tshort_contigs\n")

        for rpt_path in base_dir.glob("*_Viruses_kraken_report.txt"):
            sid = rpt_path.stem.replace("_Viruses_kraken_report", "")
            kout_path = base_dir / f"{sid}_kraken2_output.txt"
            contigs_path = base_dir / sid / "contigs.fasta"

            if not rpt_path.exists() or not kout_path.exists() or not contigs_path.exists():
                print(f"Skipping {sid}: missing required files.")
                continue

            print(f"\nExtracting from sample {sid}")
            taxon_map = {}
            with open(rpt_path) as f:
                for row in csv.reader(f, delimiter="\t"):
                    if len(row) < 6 or row[3] not in allowed_ranks:
                        continue
                    tid = row[4].strip()
                    pname = row[5].strip().replace(" ", "_").replace("\\", "_").replace("/", "_")
                    if any(tok in pname.lower() for tok in ("virus", "virinae", "viridae")):
                        taxon_map[tid] = pname
            if not taxon_map:
                continue

            contig_dict = {}
            with open(kout_path) as f:
                for row in csv.reader(f, delimiter="\t"):
                    if len(row) < 3:
                        continue
                    cid, info = row[1].strip(), row[2]
                    for tid in taxon_map:
                        if f"taxid {tid}" in info:
                            contig_dict.setdefault(tid, []).append(cid)
            if not contig_dict:
                continue

            for tid, clist in contig_dict.items():
                pname = taxon_map[tid]
                outfa = contigs_path.parent / f"{sid}_{pname}.fasta"
                long_c = short_c = 0
                with open(outfa, "w") as out:
                    for rec in SeqIO.parse(str(contigs_path), "fasta"):
                        if rec.id in clist:
                            if len(rec.seq) > 500:
                                SeqIO.write(rec, out, "fasta")
                                long_c += 1
                            else:
                                short_c += 1

                if long_c == 0:
                    outfa.unlink()
                    print(f"  Deleted {outfa} ({short_c} ≤500 bp)")
                else:
                    print(f"  Wrote {outfa}: {long_c} >500 bp, {short_c} ≤500 bp")
                summary_f.write(f"{sid}\t{tid}\t{pname}\t{long_c}\t{short_c}\n")

    print(f"\nContig extraction done; summary → {summary_filename}")

def merge_and_rename_contigs(base_contigs_dir, merged_filename="merged_contigs_renamed.fasta"):
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
