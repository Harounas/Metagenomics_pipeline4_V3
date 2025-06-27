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
    """Run RagTag to scaffold *query* contigs against a *reference* genome.
    Removes *output_dir* if it already exists, then executes `ragtag.py scaffold`.
    Returns True on success, False otherwise.
    """
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
    """Write only the first FASTA record from *input_fasta* to *output_fasta* using seqtk/awk."""
    cmd = (
        f"seqtk seq -A {input_fasta} | "
        "awk '/^>/ {if (++c>1) exit} {print}' > " + output_fasta
    )
    try:
        subprocess.run(cmd, shell=True, check=True)
        return True
    except subprocess.CalledProcessError:
        return False


def fetch_reference_from_ncbi(acc: str, out_fasta: str, email: str = "you@example.com") -> str:
    """Download a reference genome FASTA from NCBI (nucleotide) using an accession."""
    Entrez.email = email
    with Entrez.efetch(db="nuccore", id=acc, rettype="fasta", retmode="text") as handle:
        fasta_text = handle.read()

    with open(out_fasta, "w") as fh:
        fh.write(fasta_text)

    return out_fasta


def select_contigs(contigs_fasta: str, contig_ids: Set[str], out_fasta: str) -> int:
    """Write contigs whose IDs are in *contig_ids* to *out_fasta*. Returns number written."""
    records = [rec for rec in SeqIO.parse(contigs_fasta, "fasta") if rec.id in contig_ids]
    SeqIO.write(records, out_fasta, "fasta")
    return len(records)


def load_virus_reference_map(map_tsv: str) -> dict:
    """Return dict mapping virus name -> accession from a 2‑column TSV."""
    df = pd.read_csv(map_tsv, sep="\t")
    return dict(zip(df["virus"], df["accession"]))


def valid_base_count(seq: str) -> int:
    """Count A/C/G/T bases only (case‑insensitive)."""
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
    virus_ref_map_tsv: Optional[str] = None,
    reference_acc: Optional[str] = None,
    reference_fasta: Optional[str] = None,
    threads: int = 8,
    tmp_dir: str = "tmp_scaff",
    email: str = "you@example.com",
) -> tuple[str, str]:
    """Scaffold viral contigs and emit a summary TSV.

    Returns:
        (path_to_final_fasta, path_to_summary_tsv)
    """
    os.makedirs(tmp_dir, exist_ok=True)
    os.makedirs(output_root, exist_ok=True)

    # 1 ─────────── Filter TSV for sample & virus ──────────────────────────────
    df = pd.read_csv(tsv_path, sep="\t")
    rows = df[(df["Sample_ID"] == sample_id) & (df["virus"] == virus_name)]
    if rows.empty:
        raise ValueError(f"No rows found for sample '{sample_id}' & virus '{virus_name}'.")

    contig_ids: Set[str] = set(rows["query_id"].unique())

    # 2 ─────────── Extract contigs from per‑sample assembly ───────────────────
    sample_contigs_fasta = Path(contigs_root) / sample_id / "contigs.fasta"
    if not sample_contigs_fasta.exists():
        raise FileNotFoundError(f"Contigs FASTA not found: {sample_contigs_fasta}")

    selected_contigs_fasta = Path(tmp_dir) / f"{sample_id}_{virus_name.replace(' ', '_')}_selected.fasta"
    n_selected = select_contigs(str(sample_contigs_fasta), contig_ids, str(selected_contigs_fasta))
    if n_selected == 0:
        raise ValueError("No contigs extracted; aborting.")

    # 3 ─────────── Determine reference accession (mapping file or arg) ────────
    if not reference_fasta and not reference_acc:
        if virus_ref_map_tsv is None:
            raise ValueError("No reference provided and no virus_ref_map_tsv supplied.")
        virus_map = load_virus_reference_map(virus_ref_map_tsv)
        reference_acc = virus_map.get(virus_name)
        if reference_acc is None:
            raise ValueError(f"No accession found for virus '{virus_name}' in map file.")

    # 4 ─────────── Fetch reference if needed ──────────────────────────────────
    if reference_fasta:
        reference_path = reference_fasta
    else:
        reference_path = Path(tmp_dir) / f"{reference_acc}.fasta"
        if not reference_path.exists():
            fetch_reference_from_ncbi(reference_acc, str(reference_path), email=email)

    # 5 ─────────── RagTag scaffold ────────────────────────────────────────────
    ragtag_out = Path(tmp_dir) / f"ragtag_{sample_id}_{virus_name.replace(' ', '_')}"
    success = run_ragtag_scaffold(
        str(reference_path), str(selected_contigs_fasta), str(ragtag_out), threads=threads
    )
    if not success:
        raise RuntimeError("RagTag scaffold failed.")

    ragtag_scaffold_fasta = ragtag_out / "ragtag.scaffold.fasta"
    if not ragtag_scaffold_fasta.exists():
        raise FileNotFoundError("ragtag.scaffold.fasta not produced.")

    # 6 ─────────── Rename first scaffold & write final FASTA ──────────────────
    final_fasta = Path(output_root) / f"{sample_id}_{virus_name.replace(' ', '_')}.fasta"

    tmp_first = ragtag_out / "first_scaffold.fasta"
    extract_first_contig(str(ragtag_scaffold_fasta), str(tmp_first))

    records = list(SeqIO.parse(tmp_first, "fasta"))
    records[0].id = f"{sample_id}_{virus_name.replace(' ', '_')}"
    records[0].description = ""
    SeqIO.write(records, str(final_fasta), "fasta")

    # 7 ─────────── Write contig‑length summary TSV ────────────────────────────
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
