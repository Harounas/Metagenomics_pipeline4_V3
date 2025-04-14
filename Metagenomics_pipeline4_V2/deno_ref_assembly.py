import os
import subprocess
from Bio import SeqIO
import pandas as pd
import logging

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("pipeline.log"),
        logging.StreamHandler()
    ]
)

def download_and_index_reference(tax, scientific_name, tax_dir):
    """
    Download and index the reference genome for a given taxonomic ID.
    """
    fasta_file = os.path.join(tax_dir, f"{scientific_name}.fasta")
    command = f'esearch -db nucleotide -query "txid{tax}[Organism]" | efilter -source refseq | efetch -format fasta > {fasta_file}'
    try:
        logging.info(f"Downloading reference genome for {scientific_name} (taxid {tax})...")
        subprocess.run(command, shell=True, check=True)
        subprocess.run(f"bwa index {fasta_file}", shell=True, check=True)
        logging.info(f"Reference genome downloaded and indexed: {fasta_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error preparing reference for {tax}: {e}")
        return None
    return fasta_file

def deno_ref_based(sample, sample_r1, sample_r2, output_dir):
    """
    Run de novo assembly using MetaSPAdes.
    """
    contigs_file = os.path.join(f"{output_dir}/{sample}_denovo", "contigs.fasta")
    if os.path.exists(contigs_file):
        logging.info(f"De novo assembly already exists for {sample}. Skipping.")
        return contigs_file

    command = f"metaspades.py -1 {sample_r1} -2 {sample_r2} -o {output_dir}/{sample}_denovo -t 32"
    try:
        logging.info(f"Running de novo assembly for {sample}: {command}")
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in de novo assembly for {sample}: {e}")
        return None

    if not os.path.exists(contigs_file):
        logging.error(f"De novo assembly failed: contigs not generated for {sample}")
        return None

    logging.info(f"De novo assembly completed for {sample}: {contigs_file}")
    return contigs_file
