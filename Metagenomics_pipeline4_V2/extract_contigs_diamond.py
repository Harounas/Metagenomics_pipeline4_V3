#!/usr/bin/env python3
"""
extract_contigs_diamond.py

This script extracts viral contigs from Kraken2 results, merges & renames them,
runs geNomad, clusters contigs, extracts long contigs from the clustered FASTA,
and performs Diamond BLASTX with Entrez annotation.
"""

from pathlib import Path
import subprocess
import pandas as pd
from Bio import SeqIO, Entrez
import sys

# Set your email for NCBI Entrez usage
Entrez.email = "harounasoum17@gmail.com"


import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def extract_long_contigs_kraken(base_contigs_dir, output_tsv="long_contigs_summary.tsv") -> Path:
    """
    Extract viral contigs with length ≥500 bp from Kraken2 output.
    """
    from Bio import SeqIO

    allowed_ranks = {"F", "F1", "F2", "G", "G1", "G2", "S", "S1", "S2"}
    base_dir = Path(base_contigs_dir)
    output_dir = base_dir / "long_contigs"
    output_dir.mkdir(exist_ok=True)

    all_long_seqs = []
    header = ["Sample_ID", "gene", "contig_length", "taxname"]

    with open(output_tsv, "w", newline="") as out_f:
        writer = csv.writer(out_f, delimiter="\t")
        writer.writerow(header)

        for report in base_dir.glob("*_Viruses_kraken_report.txt"):
            sample_id = report.stem.replace("_Viruses_kraken_report", "")
            kout = base_dir / f"{sample_id}_kraken2_output.txt"
            contigs_fa = base_dir / sample_id / "contigs.fasta"

            if not (report.exists() and kout.exists() and contigs_fa.exists()):
                print(f"Skipping {sample_id}: missing files.")
                continue

            # 1. Build taxid → virus name map
            taxon_map = {}
            with open(report) as rf:
                reader = csv.reader(rf, delimiter="\t")
                for row in reader:
                    if len(row) < 6 or row[3].strip() not in allowed_ranks:
                        continue
                    taxid = row[4].strip()
                    name = row[5].strip().replace(" ", "_")
                    if any(tok in name.lower() for tok in ("virus", "virinae", "viridae")):
                        taxon_map[taxid] = name

            if not taxon_map:
                continue

            # 2. Map taxid → list of contig IDs
            contig_hits = {}
            with open(kout) as kf:
                reader = csv.reader(kf, delimiter="\t")
                for row in reader:
                    if len(row) < 3:
                        continue
                    cid, info = row[1].strip(), row[2]
                    for tid in taxon_map:
                        if f"taxid {tid}" in info:
                            contig_hits.setdefault(tid, []).append(cid)

            if not contig_hits:
                continue

            contig_dict = {rec.id: rec for rec in SeqIO.parse(str(contigs_fa), "fasta")}

            for taxid, contig_ids in contig_hits.items():
                taxname = taxon_map[taxid]
                for cid in contig_ids:
                    if cid in contig_dict:
                        rec = contig_dict[cid]
                        if len(rec.seq) >= 400:
                            new_id = f"{sample_id}_{cid}"
                            new_rec = SeqRecord(rec.seq, id=new_id, description="")
                            all_long_seqs.append(new_rec)
                            writer.writerow([sample_id, cid, len(rec.seq), taxname])

    output_fasta = output_dir / "all_long_contigs_renamed.fasta"
    with open(output_fasta, "w") as fasta_out:
        SeqIO.write(all_long_seqs, fasta_out, "fasta")

    print(f"✅ Merged long Kraken contigs saved to: {output_fasta}")
    return output_fasta

# Example usage
#extract_long_contigs(".", output_tsv="long_contigs_summary.tsv")


def extract_short_contigs_kraken(base_contigs_dir, output_tsv="short_contigs_summary.tsv"):
    """
    For each sample and each virus taxon:
      - find all contigs assigned by Kraken2
      - compute each contig’s length from contigs.fasta
      - if the largest contig length < 500 bp, emit all contigs for that virus
        into a TSV [sample_id, contig_id, contig_length, virus_name]
    Expects files under base_contigs_dir:
      - {sample_id}_Viruses_kraken_report.txt
      - {sample_id}_kraken2_output.txt
      - {sample_id}/contigs.fasta
    """
    allowed_ranks = {"F", "F1", "F2", "G", "G1", "G2", "S", "S1", "S2"}
    base_dir = Path(base_contigs_dir)
    header = ["Sample_ID", "gene", "contig_length", "taxname"]
    
    with open(output_tsv, "w", newline="") as out_f:
        writer = csv.writer(out_f, delimiter="\t")
        writer.writerow(header)

        # iterate through each sample’s virus-only Kraken report
        for report in base_dir.glob("*_Viruses_kraken_report.txt"):
            sample_id = report.stem.replace("_Viruses_kraken_report", "")
            kout = base_dir / f"{sample_id}_kraken2_output.txt"
            contigs_fa = base_dir / sample_id / "contigs.fasta"

            # skip if any required file is missing
            if not (report.exists() and kout.exists() and contigs_fa.exists()):
                print(f"Skipping {sample_id}: missing files.")
                continue

            # 1) build taxon_id → virus_name map
            taxon_map = {}
            with open(report) as rf:
                reader = csv.reader(rf, delimiter="\t")
                for row in reader:
                    if len(row) < 6 or row[3] not in allowed_ranks:
                        continue
                    taxid = row[4].strip()
                    name = row[5].strip().replace(" ", "_")
                    if any(tok in name.lower() for tok in ("virus","virinae","viridae")):
                        taxon_map[taxid] = name

            if not taxon_map:
                continue

            # 2) build taxon_id → list of contig IDs from Kraken2 output
            contig_hits = {}
            with open(kout) as kf:
                reader = csv.reader(kf, delimiter="\t")
                for row in reader:
                    if len(row) < 3:
                        continue
                    cid, info = row[1].strip(), row[2]
                    for tid in taxon_map:
                        if f"taxid {tid}" in info:
                            contig_hits.setdefault(tid, []).append(cid)

            if not contig_hits:
                continue

            # 3) parse contigs.fasta once into a dict of lengths
            length_map = {rec.id: len(rec.seq) 
                          for rec in SeqIO.parse(str(contigs_fa), "fasta")}


            # 4) for each virus, check max length < 500
            for tid, cids in contig_hits.items():
                lengths = [length_map.get(cid, 0) for cid in cids]
                if not lengths:
                    continue
                if max(lengths) < 500:
                    vname = taxon_map[tid]
                    for cid, clen in zip(cids, lengths):
                        writer.writerow([sample_id, cid, clen, vname])
                    print(f"{sample_id} • {vname}: all contigs ≤ {max(lengths)} bp")

    print(f"\nDone! Summary written to {output_tsv}")

#extract_short_contigs(".")


def extract_and_merge_contigs_genomad(base_contigs_dir: str,
                              output_fasta: str = "merged_contigs_genomad.fasta",
                              min_length: int = 200) -> None:

    #Extracts contigs > min_length bp from each sample's contigs.fasta
    #under base_contigs_dir and merges them into a single FASTA.

    base_dir = Path(base_contigs_dir)
    merged_records = []

    for contigs_file in base_dir.glob("*/contigs.fasta"):
        sample_id = contigs_file.parent.name
        print(f"Processing sample: {sample_id}")
        for rec in SeqIO.parse(contigs_file, "fasta"):
            if len(rec.seq) > min_length:
                rec.id = f"{sample_id}|{rec.id}"
                rec.description = ""
                merged_records.append(rec)

    if merged_records:
        SeqIO.write(merged_records, output_fasta, "fasta")
        print(f"✅ Merged {len(merged_records)} contigs >{min_length} bp into {output_fasta}")
    else:
        print(f"⚠️ No contigs found >{min_length} bp in any sample.")





def run_genomad(input_fasta: str, output_dir: str, genomad_db: str,
                min_score: float = 0.5, threads: int = 8) -> Path:
    """
    Runs geNomad end-to-end pipeline to identify viral contigs.
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    genomad_cmd = [
    "genomad", "end-to-end",
    input_fasta,
    str(out_dir),
    genomad_db,
    "--min-score", str(min_score),
    "--threads", str(threads),
        "--restart"
]
    print("🚀 Running geNomad:\n", " ".join(genomad_cmd))
    subprocess.run(genomad_cmd, check=True)

    virus_fasta = out_dir / "merged_contigs_genomad_summary" / "merged_contigs_genomad_virus.fna"
    if not virus_fasta.exists():
        raise FileNotFoundError(f"Viral contigs FASTA not found: {virus_fasta}")

    return virus_fasta


def filter_and_merge(fasta_paths, min_length, output_path):
    seen_ids = set()
    records_to_write = []

    for fasta in fasta_paths:
        for rec in SeqIO.parse(str(fasta), "fasta"):
            if len(rec.seq) >= min_length:
                if rec.id not in seen_ids:
                    seen_ids.add(rec.id)
                    records_to_write.append(rec)

    SeqIO.write(records_to_write, output_path, "fasta")
    print(f"✅ Wrote {len(records_to_write)} contigs ≥{min_length} bp to {output_path}")

def cluster_contigs(virus_fasta: Path, output_dir: str, final_output: str = "clustered_contigs.fasta",
                    identity: float = 0.95, word_size: int = 10,
                    mem_mb: int = 16000, threads: int = 8) -> Path:
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    clustered_fasta = out_dir / final_output
    cdhit_cmd = [
        "cd-hit-est",
        "-i", str(virus_fasta),
        "-o", str(clustered_fasta),
        "-c", str(identity),
        "-n", str(word_size),
        "-d", "0",
        "-M", str(mem_mb),
        "-T", str(threads)
    ]
    print("🚀 Running cd-hit-est:\n", " ".join(cdhit_cmd))
    subprocess.run(cdhit_cmd, check=True)

    return clustered_fasta
"""
def run_genomad_and_cluster(input_fasta: str,
                            output_dir: str,
                            genomad_db: str,
                            final_output: str = "clustered_contigs.fasta",
                            min_score: float = 0.5,
                            identity: float = 0.95,
                            word_size: int = 10,
                            mem_mb: int = 16000,
                            threads: int = 8) -> None:
    
    #Runs geNomad to identify viral contigs, then clusters those contigs
    #using CD-HIT-EST. Produces a clustered FASTA at output_dir/final_output.
    
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # geNomad end-to-end
    genomad_cmd = [
        "genomad", "end-to-end",
        input_fasta,
        str(out_dir),
        "--db-dir", genomad_db,
        "--min-score", str(min_score),
        "--threads", str(threads)
    ]
    print("🚀 Running geNomad:\n", " ".join(genomad_cmd))
    subprocess.run(genomad_cmd, check=True)

    # Locate filtered viral contigs
    virus_fasta = out_dir / "merged_contigs_filtered_summary" / "merged_contigs_virus.fna"
    if not virus_fasta.exists():
        raise FileNotFoundError(f"Viral contigs FASTA not found: {virus_fasta}")

    # Cluster with CD-HIT-EST
    clustered_fasta = out_dir / final_output
    cdhit_cmd = [
        "cd-hit-est",
        "-i", str(virus_fasta),
        "-o", str(clustered_fasta),
        "-c", str(identity),
        "-n", str(word_size),
        "-d", "0",
        "-M", str(mem_mb),
        "-T", str(threads)
    ]
    print("🚀 Running cd-hit-est:\n", " ".join(cdhit_cmd))
    subprocess.run(cdhit_cmd, check=True)

    print(f"\n✅ Pipeline complete! Final clustered FASTA: {clustered_fasta}")

"""

def extract_long_contigs(input_fasta: str,
                          output_fasta: str,
                          min_length: int = 400) -> None:
    """
    Extracts contigs longer than min_length from the clustered contigs FASTA
    produced by run_genomad_and_cluster. Writes output to output_fasta.
    """
    records = []
    for rec in SeqIO.parse(input_fasta, "fasta"):
        if len(rec.seq) > min_length:
            records.append(rec)

    if records:
        SeqIO.write(records, output_fasta, "fasta")
        print(f"✅ Extraction complete: {len(records)} contigs >{min_length} bp written to {output_fasta}")
    else:
        print(f"⚠️ No contigs >{min_length} bp in {input_fasta}")



def run_diamond(diamond_db: str,
                query_file: str,
                output_file: str = "results.m8",
                threads: int = 8) -> None:
    """
    Runs Diamond BLASTX of query_file against diamond_db.
    """
    cmd = [
        "diamond", "blastx",
        "--query", query_file,
        "--db", diamond_db,
        "--out", output_file,
        "--threads", str(threads),
        "--outfmt", "6"  # tabular
    ]
    print("\n🚀 Executing Diamond BLASTX:\n", " ".join(cmd))
    subprocess.run(cmd, check=True)
    print(f"✅ Diamond results written to {output_file}")

def process_diafmond_results(results_file: str = "results.m8",
                            out_csv: str = "extracted_virus.csv",
                            sorted_csv: str = "extracted_virus_sorted.csv") -> None:
    """
    Parses Diamond BLASTX output, annotates via Entrez, computes a custom score,
    and writes both raw and sorted CSV results.
    """
    df = pd.read_csv(results_file, sep="\t", header=None, names=[
        'query_id','subject_id','pident','aln_len',
        'mismatches','gaps','qstart','qend','sstart','send',
        'evalue','bitscore'
    ])

    def fetch_name(acc: str) -> str:
        try:
            handle = Entrez.efetch(db="protein", id=acc, rettype="gb", retmode="text")
            txt = handle.read()
            handle.close()
            for line in txt.splitlines():
                if line.startswith("  ORGANISM"):
                    return line.split("ORGANISM")[-1].strip()
        except Exception:
            return "Unknown"
        return "Unknown"

    print("🔍 Annotating subject IDs via Entrez...")
    df['Virus Name'] = df['subject_id'].apply(fetch_name)
    df.to_csv(out_csv, index=False)

    df['Custom_Score'] = (df.pident * df.aln_len) / (1 + df.mismatches + df.gaps)
    df = df.sort_values('Custom_Score', ascending=False)
    df.to_csv(sorted_csv, index=False)

    print(f"✅ Annotated Diamond output → {out_csv}")
    print(f"✅ Sorted results → {sorted_csv}")
