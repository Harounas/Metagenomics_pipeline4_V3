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
import os
import re

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
                        if len(rec.seq) >= 200:
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
"""

def extract_short_contigs_kraken(base_contigs_dir, output_tsv="short_contigs_summary.tsv"):
    """
    #For each sample and each virus taxon:
      #- find all contigs assigned by Kraken2
      #- compute each contig’s length from contigs.fasta
      #- if the largest contig length < 500 bp, emit all contigs for that virus
     #   into a TSV [sample_id, contig_id, contig_length, virus_name]
    #Expects files under base_contigs_dir:
     # - {sample_id}_Viruses_kraken_report.txt
    #  - {sample_id}_kraken2_output.txt
   #   - {sample_id}/contigs.fasta
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
"""

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



def save_genomad_short_virus_contigs(
    input_tsv: str,
    output_tsv: str = "short_contig_virus_table.tsv",
    max_virus_len: int = 500,   # virus kept only if its longest contig < 500 bp
    min_contig_len: int = 200,  # discard any contig shorter than this
    suffix_re: str = r"(.+?)_\d+$",  # strip trailing “_1”, “_2”, …
):
    """
    From a geNomad merged-genes TSV create a 4-column table:
        Sample_ID   Sample_ID|gene   contig_length   taxname

    • A virus (sample × taxname pair) is included only if *all* its contigs
      are < `max_virus_len` bp.
    • Individual contigs shorter than `min_contig_len` are discarded.
    • Trailing “_N” suffixes are removed from the contig part of the gene.
    • Underscores in taxname are replaced with spaces.
    """

    suffix_pat      = re.compile(suffix_re)
    contig_len_pat  = re.compile(r"length_(\d+)")      # pulls “230” from “…length_230_cov_…”
    # 1️⃣  First pass – collect max-contig length per (sample, taxname)
    max_len = defaultdict(int)   # (sample_id, taxname) ➜ longest contig length

    rows_tmp = []                # will hold rows for second pass
    with open(input_tsv, newline="") as fin:
        reader = csv.DictReader(fin, delimiter="\t")
        for row in reader:
            raw_gene  = row["gene"]
            clean_gene = suffix_pat.sub(r"\1", raw_gene)

            sample_id  = clean_gene.split("|", 1)[0]
            match_len  = contig_len_pat.search(clean_gene)
            if not match_len:
                continue                             # skip if length not found in ID
            contig_len = int(match_len.group(1))

            taxname = row.get("taxname", "").replace("_", " ")
            key = (sample_id, taxname)
            if contig_len > max_len[key]:
                max_len[key] = contig_len

            rows_tmp.append((sample_id, clean_gene, contig_len, taxname))

    # 2️⃣  Second pass – write rows that satisfy both filters
    with open(output_tsv, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(["Sample_ID", "Sample_ID|gene", "contig_length", "taxname"])

        kept = 0
        for sample_id, clean_gene, contig_len, taxname in rows_tmp:
            if contig_len < min_contig_len:
                continue                                            # too short
            if max_len[(sample_id, taxname)] >= max_virus_len:
                continue                                            # virus too long
            writer.writerow([sample_id, clean_gene, contig_len, taxname])
            kept += 1

    print(f"✓ Wrote {kept:,} contig rows to {output_tsv}")





# -----------------------------------------------------------
# 1) geNomad gene table → short-virus contigs (4-column TSV)
# -----------------------------------------------------------

def save_genomad_short_virus_contigs(
    input_tsv: str,
    output_tsv: str = "short_contig_virus_table.tsv",
    output_dir: str = ".",                 # NEW
    max_virus_len: int = 500,              # virus kept only if its longest contig < 500 bp
    min_contig_len: int = 200,             # discard any contig shorter than this
    suffix_re: str = r"(.+?)_\d+$"         # strip trailing “_1”, “_2”, …
):
    """
    Writes <output_dir>/<output_tsv> with columns:
        Sample_ID   Sample_ID|gene   contig_length   taxname
    Filters:
      • keep a virus (sample × taxname) only if *all* its contigs < max_virus_len
      • keep individual contigs only if length ≥ min_contig_len
    """
    suffix_pat      = re.compile(suffix_re)
    contig_len_pat  = re.compile(r"length_(\d+)")

    # ── 1️⃣ first pass: find longest contig per (sample, taxname)
    max_len = defaultdict(int)          # (sample, taxname) → longest length
    rows_tmp = []                       # cache rows for second pass

    with open(input_tsv, newline="") as fin:
        reader = csv.DictReader(fin, delimiter="\t")
        for row in reader:
            raw_gene   = row["gene"]
            clean_gene = suffix_pat.sub(r"\1", raw_gene)
            sample_id  = clean_gene.split("|", 1)[0]

            m_len = contig_len_pat.search(clean_gene)
            if not m_len:
                continue
            contig_len = int(m_len.group(1))

            taxname = row.get("taxname", "").replace("_", " ")
            key = (sample_id, taxname)
            max_len[key] = max(max_len[key], contig_len)

            rows_tmp.append((sample_id, clean_gene, contig_len, taxname))

    # ── 2️⃣ second pass: write filtered rows
    output_path = Path(output_dir) / output_tsv
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", newline="") as fout:
        w = csv.writer(fout, delimiter="\t")
        w.writerow(["Sample_ID", "Sample_ID|gene", "contig_length", "taxname"])

        kept = 0
        for sample_id, clean_gene, contig_len, taxname in rows_tmp:
            if contig_len < min_contig_len:
                continue
            if max_len[(sample_id, taxname)] >= max_virus_len:
                continue
            w.writerow([sample_id, clean_gene, contig_len, taxname])
            kept += 1

    print(f"✓ Wrote {kept:,} contig rows to {output_path}")



# -----------------------------------------------------------
# 2) Kraken-based filter of short-virus contigs
# -----------------------------------------------------------

def extract_short_contigs_kraken(
    base_contigs_dir,
    output_tsv="short_contigs_summary.tsv",
    output_dir: str = ".",                # NEW
    max_virus_len=500,                    # virus kept only if its longest contig < 500 bp
    min_contig_len=200,                   # discard any contig shorter than this
):
    """
    Writes <output_dir>/<output_tsv> with columns:
        Sample_ID   Sample_ID|gene   contig_length   taxname
    """
    allowed_ranks = {"S", "S1", "S2"}
    base_dir = Path(base_contigs_dir)
    output_path = Path(output_dir) / output_tsv
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", newline="") as out_f:
        writer = csv.writer(out_f, delimiter="\t")
        writer.writerow(["Sample_ID", "Sample_ID|gene", "contig_length", "taxname"])

        for report in base_dir.glob("*_Viruses_kraken_report.txt"):
            sample_id  = report.stem.replace("_Viruses_kraken_report", "")
            kout       = base_dir / f"{sample_id}_kraken2_output.txt"
            contigs_fa = base_dir / sample_id / "contigs.fasta"
            if not (kout.exists() and contigs_fa.exists()):
                print(f"⚠  Skipping {sample_id}: missing files.")
                continue

            # 1️⃣ taxid → virus name
            taxon_map = {}
            with open(report) as rf:
                for row in csv.reader(rf, delimiter="\t"):
                    if len(row) < 6 or row[3] not in allowed_ranks:
                        continue
                    taxid = row[4].strip()
                    name  = row[5]
                    if any(tok in name.lower() for tok in ("virus", "virinae", "viridae")):
                        taxon_map[taxid] = name
            if not taxon_map:
                continue

            # 2️⃣ taxid → contig IDs
            contig_hits = {}
            with open(kout) as kf:
                for row in csv.reader(kf, delimiter="\t"):
                    if len(row) < 3:
                        continue
                    cid, info = row[1].strip(), row[2]
                    for tid in taxon_map:
                        if f"taxid {tid}" in info:
                            contig_hits.setdefault(tid, []).append(cid)
            if not contig_hits:
                continue

            # 3️⃣ contig length lookup
            length_map = {rec.id: len(rec.seq) for rec in SeqIO.parse(contigs_fa, "fasta")}

            # 4️⃣ filter & write
            for tid, cids in contig_hits.items():
                lengths = [length_map.get(cid, 0) for cid in cids]
                if not lengths or max(lengths) >= max_virus_len:
                    continue

                vname = taxon_map[tid]
                for cid, clen in zip(cids, lengths):
                    if clen >= min_contig_len:
                        writer.writerow([sample_id, f"{sample_id}|{cid}", clen, vname])

    print(f"✓ Summary written to {output_path}")


def merge_short_contigs(
    tsv_a: str,
    tsv_b: str,
    output_tsv: str = "combined_no_dupes.tsv",
    output_dir: str = "."
):
    """
    Combine two TSV tables and remove duplicate (Sample_ID|gene, taxname) pairs.

    Parameters:
    -----------
    tsv_a : str
        Path to the first TSV file.
    tsv_b : str
        Path to the second TSV file.
    output_tsv : str
        Name of the output TSV file.
    output_dir : str
        Directory where output TSV will be written.

    Returns:
    --------
    Writes the combined TSV file to: {output_dir}/{output_tsv}
    """
    output_path = Path(output_dir) / output_tsv
    output_path.parent.mkdir(parents=True, exist_ok=True)  # ensure directory exists

    uniq = set()  # Track (Sample_ID|gene, taxname) pairs
    written = 0

    with open(output_path, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(["Sample_ID", "Sample_ID|gene", "contig_length", "taxname"])

        for src in (tsv_a, tsv_b):
            with open(src, newline="") as fin:
                reader = csv.reader(fin, delimiter="\t")
                header = next(reader)  # skip header
                for row in reader:
                    contig_id = row[1]
                    taxname   = row[3]
                    key = (contig_id, taxname)
                    if key in uniq:
                        continue
                    uniq.add(key)
                    writer.writerow(row)
                    written += 1

    print(f"✓ Merged:\n   {Path(tsv_a)}\n   {Path(tsv_b)}")
    print(f"→ Wrote {written:,} unique rows to {output_path}")


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





def process_virus_contigs(fasta_file, diamond_results_file, output_dir):
    """
    Processes virus annotations and DIAMOND BLASTX results to output annotated results per contig.

    Args:
        fasta_file (str): Path to the input FASTA file (e.g. "MNT/nr.fasta")
        diamond_results_file (str): DIAMOND results .m8 format file path
        output_dir (str): Directory to write outputs
    """
    os.makedirs(output_dir, exist_ok=True)
    virus_list_path = os.path.join(output_dir, "virus_list.tsv")

    # Step 1: Generate virus list only if not already present
    if not os.path.exists(virus_list_path):
        print("Generating virus list from FASTA...")
        pattern = re.compile(r'^>(\S+).*?\[([^\]]*virus[^\]]*)\]', re.IGNORECASE)
        with open(fasta_file, 'r') as infile, open(virus_list_path, 'w') as outfile:
            outfile.write("accession_number\tvirus_name\n")
            for line in infile:
                match = pattern.match(line)
                if match:
                    accession = match.group(1)
                    virus_name = match.group(2)
                    outfile.write(f"{accession}\t{virus_name}\n")
        print(f"Virus entries extracted to: {virus_list_path}")
    else:
        print(f"Virus list already exists at: {virus_list_path}, skipping generation.")

    # Step 2: Process DIAMOND results and extract best hits per query
    df = pd.read_csv(diamond_results_file, sep='\t', header=None)
    best_hits = df.loc[df.groupby(0)[11].idxmax()]

    # Step 3: Load virus list and map subject IDs to virus names
    virus_df = pd.read_csv(virus_list_path, sep='\t')
    virus_map = dict(zip(virus_df['accession_number'], virus_df['virus_name']))
    best_hits.columns = [
        "query_id", "subject_id", "pident", "aln_len", "mismatches", "gaps",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    ]
    best_hits['virus'] = best_hits['subject_id'].map(virus_map)
    best_hits.drop(columns=["subject_id"], inplace=True)

    # Step 4: Parse metadata from query ID
    best_hits['Sample_ID'] = best_hits['query_id'].str.split('|').str[0]
    contig_part = best_hits['query_id'].str.split('|').str[-1]
    best_hits['contigs_len'] = pd.to_numeric(contig_part.str.split('_').str[3], errors='coerce')

    # Step 5: Filter contigs >= 500 and reorder columns
    filtered_df = best_hits[best_hits['contigs_len'] >= 500]
    col_order = ['Sample_ID', 'contigs_len'] + [
        "query_id", "pident", "aln_len", "mismatches", "gaps", "qstart", "qend",
        "sstart", "send", "evalue", "bitscore", "virus"
    ]
    filtered_df = filtered_df[col_order]

    # Step 6: Write the result
    output_path = os.path.join(output_dir, 'diamond_results_contig_with_sampleid.tsv')
    filtered_df.to_csv(output_path, sep='\t', index=False)
    print(f"Final annotated results saved to: {output_path}")

    return output_path
def process_diamond_results(results_file: str = "results.m8",
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
