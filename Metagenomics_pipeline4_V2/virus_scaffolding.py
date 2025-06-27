import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional, Set

import pandas as pd
from Bio import Entrez, SeqIO

# -----------------------------
# Helper utilities
# -----------------------------

def run_ragtag_scaffold(reference: str, query: str, output_dir: str, threads: int = 8) -> bool:
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    try:
        subprocess.run(
            ["ragtag.py", "scaffold", "-t", str(threads), reference, query, "-o", output_dir],
            check=True,
        )
        return True
    except subprocess.CalledProcessError:
        return False

def extract_first_contig(input_fasta: str, output_fasta: str) -> bool:
    cmd = (
        f"seqtk seq -A {input_fasta} | "
        "awk '/^>/ {if (++c>1) exit} {print}' > " + output_fasta
    )
    try:
        subprocess.run(cmd, shell=True, check=True)
        return True
    except subprocess.CalledProcessError:
        return False

def fetch_reference_from_ncbi_by_name(virus_name: str, out_fasta: str, email: str = "you@example.com") -> str:
    Entrez.email = email
    search = Entrez.esearch(db="nuccore", term=f"{virus_name}[Organism] AND complete genome", retmax=1)
    record = Entrez.read(search)
    ids = record["IdList"]
    if not ids:
        raise ValueError(f"No NCBI entry found for virus name: {virus_name}")
    acc = ids[0]
    with Entrez.efetch(db="nuccore", id=acc, rettype="fasta", retmode="text") as handle:
        fasta_text = handle.read()
    with open(out_fasta, "w") as fh:
        fh.write(fasta_text)
    return out_fasta

def select_contigs(contigs_fasta: str, contig_ids: Set[str], out_fasta: str) -> int:
    records = [rec for rec in SeqIO.parse(contigs_fasta, "fasta") if rec.id in contig_ids]
    SeqIO.write(records, out_fasta, "fasta")
    return len(records)

def valid_base_count(seq: str) -> int:
    return sum(1 for b in seq if b.upper() in {"A", "C", "G", "T"})

# -----------------------------
# Core workflow
# -----------------------------

def scaffold_virus_contigs(
    tsv_path: str,
    sample_id: str,
    virus_name: str,
    contigs_root: str,
    output_root: str,
    reference_fasta: Optional[str] = None,
    threads: int = 8,
    tmp_dir: str = "tmp_scaff",
    email: str = "you@example.com",
) -> tuple[str, str]:
    os.makedirs(tmp_dir, exist_ok=True)
    os.makedirs(output_root, exist_ok=True)

    df = pd.read_csv(tsv_path, sep="\t")
    rows = df[(df["Sample_ID"] == sample_id) & (df["virus"] == virus_name)]
    if rows.empty:
        raise ValueError(f"No rows found for sample '{sample_id}' & virus '{virus_name}'.")

    contig_ids: Set[str] = set(rows["query_id"].unique())

    sample_contigs_fasta = Path(contigs_root) / sample_id / "contigs.fasta"
    if not sample_contigs_fasta.exists():
        raise FileNotFoundError(f"Contigs FASTA not found: {sample_contigs_fasta}")

    selected_contigs_fasta = Path(tmp_dir) / f"{sample_id}_{virus_name.replace(' ', '_')}_selected.fasta"
    n_selected = select_contigs(str(sample_contigs_fasta), contig_ids, str(selected_contigs_fasta))
    if n_selected == 0:
        raise ValueError("No contigs extracted; aborting.")

    if reference_fasta:
        reference_path = reference_fasta
    else:
        reference_path = Path(tmp_dir) / f"{virus_name.replace(' ', '_')}.fasta"
        if not reference_path.exists():
            fetch_reference_from_ncbi_by_name(virus_name, str(reference_path), email=email)

    ragtag_out = Path(tmp_dir) / f"ragtag_{sample_id}_{virus_name.replace(' ', '_')}"
    success = run_ragtag_scaffold(str(reference_path), str(selected_contigs_fasta), str(ragtag_out), threads=threads)
    if not success:
        raise RuntimeError("RagTag scaffold failed.")

    ragtag_scaffold_fasta = ragtag_out / "ragtag.scaffold.fasta"
    if not ragtag_scaffold_fasta.exists():
        raise FileNotFoundError("ragtag.scaffold.fasta not produced.")

    final_fasta = Path(output_root) / f"{sample_id}_{virus_name.replace(' ', '_')}.fasta"

    tmp_first = ragtag_out / "first_scaffold.fasta"
    extract_first_contig(str(ragtag_scaffold_fasta), str(tmp_first))

    records = list(SeqIO.parse(tmp_first, "fasta"))
    records[0].id = f"{sample_id}_{virus_name.replace(' ', '_')}"
    records[0].description = ""
    SeqIO.write(records, str(final_fasta), "fasta")

    total_len_valid = 0
    for rec in SeqIO.parse(str(selected_contigs_fasta), "fasta"):
        total_len_valid += valid_base_count(str(rec.seq))

    summary_tsv = Path(output_root) / "virus_contig_lengths.tsv"
    header_needed = not summary_tsv.exists()
    with open(summary_tsv, "a") as fh:
        if header_needed:
            fh.write("Sample_ID\tvirus\tcontigs_len_valid\n")
        fh.write(f"{sample_id}\t{virus_name}\t{total_len_valid}\n")

    return str(final_fasta), str(summary_tsv)
