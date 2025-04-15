#!/usr/bin/env python3
"""
kraken_abundance_pipeline.py (Version 2)

This module processes Kraken2 reports, generates abundance plots, aggregates results (with metadata or sample IDs),
and supports quality control via MultiQC. It supports de novo assembly via MetaSPAdes, optional host depletion via Bowtie2,
and allows skipping preprocessing steps.
"""

import os
import glob
import pandas as pd
import logging
import subprocess
import random
from collections import defaultdict
import plotly.express as px

# Local imports
from .trimmomatic import run_trimmomatic
from .metaspades import run_spades
from .bowtie2 import run_bowtie2
from .kraken2 import run_kraken2

# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)

def process_sample(forward, reverse, base_name, bowtie2_index, kraken_db, output_dir, threads,
                   run_bowtie, use_precomputed_reports, use_assembly,
                   skip_preprocessing=False, skip_existing=False):
    try:
        kraken_report = os.path.join(output_dir, f"{base_name}_kraken_report.txt")
        output_report = os.path.join(output_dir, f"{base_name}_kraken2_output.txt")

        if use_precomputed_reports:
            if not os.path.exists(kraken_report) or not os.path.exists(output_report):
                raise FileNotFoundError(f"Precomputed Kraken2 report or output report not found for {base_name}")
            return kraken_report, output_report

        if skip_preprocessing:
            contigs_file = os.path.join(output_dir, f"{base_name}_contigs.fasta")
            if skip_existing and os.path.exists(contigs_file):
                logging.info(f"[SKIP] Contigs exist for {base_name}, skipping assembly.")
            else:
                logging.info(f"Running SPAdes for sample {base_name}")
                contigs_file = run_spades(forward, reverse, base_name, output_dir, threads)
            kraken_input = contigs_file
        else:
            trimmed_forward = os.path.join(output_dir, f"{base_name}_1_trimmed_paired.fq.gz")
            trimmed_reverse = os.path.join(output_dir, f"{base_name}_2_trimmed_paired.fq.gz")

            if not skip_existing or not (os.path.exists(trimmed_forward) and os.path.exists(trimmed_reverse)):
                logging.info(f"Running Trimmomatic for sample {base_name}")
                run_trimmomatic(forward, reverse, base_name, output_dir, threads)

            unmapped_r1, unmapped_r2 = trimmed_forward, trimmed_reverse
            if run_bowtie:
                bowtie_unmapped_r1 = os.path.join(output_dir, f"{base_name}_1_unmapped.fq.gz")
                bowtie_unmapped_r2 = os.path.join(output_dir, f"{base_name}_2_unmapped.fq.gz")
                if not skip_existing or not (os.path.exists(bowtie_unmapped_r1) and os.path.exists(bowtie_unmapped_r2)):
                    logging.info(f"Running Bowtie2 for sample {base_name}")
                    run_bowtie2(trimmed_forward, trimmed_reverse, base_name, bowtie2_index, output_dir, threads)
                unmapped_r1, unmapped_r2 = bowtie_unmapped_r1, bowtie_unmapped_r2

            kraken_input = run_spades(unmapped_r1, unmapped_r2, base_name, output_dir, threads) if use_assembly else unmapped_r1

        if not skip_existing or not os.path.exists(kraken_report):
            logging.info(f"Running Kraken2 for sample {base_name}")
            if use_assembly or skip_preprocessing:
                run_kraken2(kraken_input, None, base_name, kraken_db, output_dir, threads)
            else:
                run_kraken2(unmapped_r1, unmapped_r2, base_name, kraken_db, output_dir, threads)

        # Process output report as well
        if not skip_existing or not os.path.exists(output_report):
            logging.info(f"Running Output Analysis for sample {base_name}")
            process_output_report(output_report, output_dir)

        return kraken_report, output_report

    except Exception as e:
        logging.error(f"Error processing sample {base_name}: {e}")
        return None, None

def process_output_report(output_report, output_dir):
    """
    Processes output report files by splitting them into domain-specific files.
    Output files are named as: {sample}_{DomainWithoutSpaces}_kraken2_output.txt.
    """
    domain_labels = {'Viruses', 'Eukaryota', 'Bacteria', 'Archaea'}
    try:
        with open(output_report, 'r') as file:
            lines = file.readlines()

        current_domain = None
        current_rows = []
        for line in lines:
            columns = line.strip().split("\t")
            if len(columns) < 6:
                continue
            rank_code = columns[3]

            if rank_code == "D":
                if current_domain:
                    save_domain_data(current_domain, current_rows, output_dir, is_kraken2_output=True)
                current_domain = columns[5]  # Domain name (e.g., Viruses, Eukaryota)
                current_rows = [line]
            else:
                current_rows.append(line)

        # Save the last domain data
        if current_domain:
            save_domain_data(current_domain, current_rows, output_dir, is_kraken2_output=True)

    except Exception as e:
        logging.error(f"Error processing output report {output_report}: {e}")

def save_domain_data(domain, rows, output_dir, is_kraken2_output=False):
    """
    Saves domain-specific data into a file.
    If is_kraken2_output is True, saves it as a Kraken2 output file.
    """
    domain_file_name = f"{domain.replace(' ', '')}_kraken2_output.txt" if is_kraken2_output else f"{domain.replace(' ', '')}_kraken_report.txt"
    domain_file_path = os.path.join(output_dir, domain_file_name)

    with open(domain_file_path, 'w') as f:
        for row in rows:
            f.write(row)

    logging.info(f"Saved {domain} data to {domain_file_path}")

def generate_sample_ids_csv(kraken_dir):
    """
    Generates a CSV containing sample IDs extracted from Kraken report filenames.

    Parameters:
      kraken_dir (str): Directory with Kraken report files.

    Returns:
      str: Path to the generated sample_ids.csv file.
    """
    try:
        sample_ids = [fname.replace('_kraken_report.txt', '')
                      for fname in os.listdir(kraken_dir)
                      if fname.endswith('_kraken_report.txt')]
        sample_ids_df = pd.DataFrame({'Sample_ID': sample_ids})
        csv_path = os.path.join(kraken_dir, 'sample_ids.csv')
        sample_ids_df.to_csv(csv_path, index=False)
        logging.info(f"Sample IDs written to {csv_path}")
        return csv_path
    except Exception as e:
        logging.error(f"Error generating sample IDs CSV: {e}")
        return None


def process_kraken_reports(kraken_dir):
    """
    Processes Kraken2 report files by splitting them into domain-specific files.
    Output files will be named as: {sample_name}_{domain}_kraken_report.txt.
    """
    domain_labels = {'Viruses', 'Eukaryota', 'Bacteria', 'Archaea'}
    
    for file_name in os.listdir(kraken_dir):
        if file_name.endswith("_kraken_report.txt"):
            kraken_report_path = os.path.join(kraken_dir, file_name)
            sample_name = clean_sample_name(file_name, domain_labels)  # Clean sample name to remove domain labels
            df = pd.read_csv(kraken_report_path, sep="\t", header=None)

            # Define the new column names
            df.columns = ["Percentage", "Reads_Covered", "Reads_Assigned", "Rank_Code", "NCBI_TaxID", "Scientific_Name"]

            # Process by domain
            for domain in domain_labels:
                # Filter the DataFrame for each domain (we assume you have a way to distinguish domains)
                domain_df = df[df['Scientific_Name'].str.contains(domain, case=False, na=False)]
                
                if not domain_df.empty:
                    # Construct new filename with sample name and domain
                    domain_output_filename = f"{sample_name}_{domain}_kraken_report.txt"
                    domain_output_path = os.path.join(kraken_dir, domain_output_filename)
                    
                    # Save the filtered data to the new file
                    domain_df.to_csv(domain_output_path, sep="\t", index=False, header=False)
                    logging.info(f"Saved {domain} data to {domain_output_path}")

def process_output_reports(kraken_dir):
    """
    Processes Kraken2 output files by matching taxon IDs and splitting them into domain-specific files.
    """
    # Define the domain labels
    domain_labels = ['Viruses', 'Eukaryota', 'Bacteria', 'Archaea']

    # Create a dictionary to store taxon IDs for each domain
    domain_taxids = {
        'Viruses': set(),
        'Bacteria': set(),
        'Eukaryota': set(),
        'Archaea': set()
    }

    # Process Kraken report files to extract taxon IDs for each domain
    for file_name in os.listdir(kraken_dir):
        if file_name.endswith("_kraken_report.txt"):  # Process Kraken report files
            domain_name = None
            if "Viruses" in file_name:
                domain_name = "Viruses"
            elif "Bacteria" in file_name:
                domain_name = "Bacteria"
            elif "Eukaryota" in file_name:
                domain_name = "Eukaryota"
            elif "Archaea" in file_name:
                domain_name = "Archaea"

            if domain_name:
                # Read the domain Kraken report to extract taxon IDs
                domain_report_path = os.path.join(kraken_dir, file_name)
                try:
                    with open(domain_report_path, 'r') as file:
                        for line in file:
                            fields = line.strip().split("\t")
                            if len(fields) >= 5:
                                taxid = fields[4]  # Extract taxid from the report
                                domain_taxids[domain_name].add(taxid)
                except Exception as e:
                    logging.error(f"Error reading {domain_name} Kraken report: {e}")

    # Now process Kraken2 output files and match taxon IDs
    for file_name in os.listdir(kraken_dir):
        if file_name.endswith("_kraken2_output.txt"):  # Kraken2 output file
            output_report_path = os.path.join(kraken_dir, file_name)
            sample_name = clean_sample_name(file_name, domain_labels)

            # Read Kraken2 output report into DataFrame
            df = pd.read_csv(output_report_path, sep="\t", header=None)

            # Ensure correct column names
            if df.shape[1] == 6:
                df.columns = ["Classification", "Node", "TaxID", "Length", "Coverage", "Taxa"]
            elif df.shape[1] == 5:
                df.columns = ["Classification", "Node", "TaxID", "Length", "Coverage"]
            else:
                logging.error(f"Unexpected number of columns in {file_name}. Skipping this file.")
                continue

            # Process and match taxon IDs for each domain
            for domain, taxids in domain_taxids.items():
                matched_taxa_df = df[df["TaxID"].isin(taxids)]

                if not matched_taxa_df.empty:
                    domain_output_filename = f"{sample_name}_kraken2_{domain}_output_report.txt"
                    domain_output_path = os.path.join(kraken_dir, domain_output_filename)
                    matched_taxa_df.to_csv(domain_output_path, sep="\t", index=False)
                    logging.info(f"Saved matched {domain} data to {domain_output_path}")

# Helper function to clean sample name
def clean_sample_name(file_name, domain_labels):
    """
    Removes domain labels and the _report.txt suffix from a filename to obtain a clean sample name.
    """
    name = file_name.replace("_report.txt", "")
    parts = name.split("_")
    cleaned_parts = [p for p in parts if p not in domain_labels]
    return "_".join(cleaned_parts)
