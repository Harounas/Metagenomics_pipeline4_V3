#!/usr/bin/env python3
"""
kraken_abundance_pipeline.py (Version 2)

This module processes Kraken2 reports, generates abundance plots, aggregates results (with metadata or sample IDs), and supports quality control via MultiQC. It supports de novo assembly via MetaSPAdes, optional host depletion via Bowtie2, and allows skipping preprocessing steps.
"""

import os
import glob
import pandas as pd
import logging
import subprocess
import random
from collections import defaultdict
import plotly.express as px
from .fastp import run_fastp

# Local imports
#from .trimmomatic import run_trimmomatic
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
        output_report = kraken_report

        if use_precomputed_reports:
            if not os.path.exists(kraken_report):
                raise FileNotFoundError(f"Precomputed Kraken2 report not found for {base_name}")
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
            trimmed_forward = os.path.join(output_dir, f"{base_name}_trimmed_R1.fastq.gz")
            trimmed_reverse = os.path.join(output_dir, f"{base_name}_trimmed_R2.fastq.gz")

            if not skip_existing or not (os.path.exists(trimmed_forward) and os.path.exists(trimmed_reverse)):
                logging.info(f"Running Trimmomatic for sample {base_name}")
                run_fastp(forward, reverse, base_name, output_dir, threads)


            unmapped_r1, unmapped_r2 = trimmed_forward, trimmed_reverse
            if run_bowtie:
                bowtie_unmapped_r1 = os.path.join(output_dir, f"{base_name}_unmapped_1.fastq.gz")
                bowtie_unmapped_r2 = os.path.join(output_dir, f"{base_name}_unmapped_2.fastq.gz")
                if not skip_existing or not (os.path.exists(bowtie_unmapped_r1) and os.path.exists(bowtie_unmapped_r2)):
                    logging.info(f"Running Bowtie2 for sample {base_name}")
                    run_bowtie2(trimmed_forward, trimmed_reverse, base_name, bowtie2_index, output_dir, threads)
                unmapped_r1, unmapped_r2 = bowtie_unmapped_r1, bowtie_unmapped_r2

            if use_assembly:
                contigs_file = os.path.join(output_dir, f"{base_name}_contigs.fasta")
                if skip_existing and os.path.exists(contigs_file):
                    logging.info(f"[SKIP] Contigs exist for {base_name}, skipping assembly.")
                else:
                    logging.info(f"Running SPAdes for sample {base_name}")
                    contigs_file = run_spades(unmapped_r1, unmapped_r2, base_name, output_dir, threads)
                kraken_input = contigs_file
            else:
                kraken_input = unmapped_r1

        if not skip_existing or not os.path.exists(kraken_report):
            logging.info(f"Running Kraken2 for sample {base_name}")
            if use_assembly or skip_preprocessing:
                run_kraken2(kraken_input, None, base_name, kraken_db, output_dir, threads)
            else:
                run_kraken2(unmapped_r1, unmapped_r2, base_name, kraken_db, output_dir, threads)

        if os.path.exists(output_report):
            logging.info(f"Running Output Analysis for sample {base_name}")
            process_output_report(output_report, output_dir)
        else:
            logging.warning(f"Skipping output analysis: {output_report} not found")

        return kraken_report, output_report

    except Exception as e:
        logging.error(f"Error processing sample {base_name}: {e}")
        return None, None


def process_output_report(output_report, output_dir):
    """
    Processes output report files by splitting them into domain-specific files.
    Output files are named as: {DomainWithoutSpaces}_output_report.txt.
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
                    save_domain_data(current_domain, current_rows, output_dir)
                current_domain = columns[5]
                current_rows = [line]
            else:
                current_rows.append(line)

        if current_domain:
            save_domain_data(current_domain, current_rows, output_dir)

    except Exception as e:
        logging.error(f"Error processing output report {output_report}: {e}")


def save_domain_data(domain, rows, output_dir):
    domain_file_path = os.path.join(output_dir, f"{domain.replace(' ', '')}_output_report.txt")
    with open(domain_file_path, 'w') as f:
        f.writelines(rows)
    logging.info(f"Saved {domain} data to {domain_file_path}")

# --- rest of functions unchanged ---

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


def remove_first_kraken(filename):
    """
    Remove the first occurrence of 'kraken' in a filename split by underscores.
    """
    parts = filename.split('_')
    if parts.count('kraken') > 1:
        idx = parts.index('kraken')
        del parts[idx]
        return '_'.join(parts)
    return filename


def process_kraken_reports(kraken_dir):
    """
    Splits Kraken summary reports into domain-specific files.
    Output files: {sample}_{Domain}_kraken_report.txt
    """
    domain_labels = {'Viruses', 'Eukaryota', 'Bacteria', 'Archaea'}
    for file_name in os.listdir(kraken_dir):
        if file_name.endswith('_kraken_report.txt'):
            path = os.path.join(kraken_dir, file_name)
            sample = clean_sample_name(file_name, domain_labels)
            domains = extract_domains_from_kraken_report(path)
            for domain, df in domains.items():
                out_name = f"{sample}_{domain.replace(' ', '')}_kraken_report.txt"
                out_name = remove_first_kraken(out_name)
                out_path = os.path.join(kraken_dir, out_name)
                df.to_csv(out_path, sep='	', index=False, header=False)
                logging.info(f"Saved {domain} data to {out_path}")


def process_output_reports(kraken_dir):
    """
    Splits domain-specific output reports (from process_output_report step).
    """
    domain_labels = {'Viruses', 'Eukaryota', 'Bacteria', 'Archaea'}
    for file_name in os.listdir(kraken_dir):
        if file_name.endswith('_output_report.txt'):
            path = os.path.join(kraken_dir, file_name)
            process_output_report(path, kraken_dir)


def clean_sample_name(file_name, domain_labels):
    """
    Removes domain labels and suffix to get base sample name.
    """
    name = file_name.replace('_report.txt', '')
    parts = name.split('_')
    cleaned = [p for p in parts if p not in domain_labels]
    return '_'.join(cleaned)


def extract_domains_from_kraken_report(kraken_report_path):
    """
    Reads a Kraken summary report and returns a dict of DataFrames by domain.
    """
    cols = ['Percentage','Reads_Covered','Reads_Assigned','Rank_Code','TaxID','Scientific_Name']
    df = pd.read_csv(kraken_report_path, sep='	', header=None, names=cols)
    domains = {}
    cd = None
    rows = []
    for _, row in df.iterrows():
        if row['Rank_Code'] == 'D':
            if cd:
                domains[cd] = pd.DataFrame(rows)
            cd = row['Scientific_Name']
            rows = [row]
        else:
            rows.append(row)
    if cd:
        domains[cd] = pd.DataFrame(rows)
    return domains


def aggregate_kraken_results(kraken_dir, metadata_file=None, sample_id_df=None,
                             read_count=1, max_read_count=10**30):
    """
    Aggregates species-level Kraken summary across samples and merges metadata.
    """
    try:
        if metadata_file:
            meta = pd.read_csv(metadata_file)
            logging.info("Using metadata file.")
        elif sample_id_df is not None:
            meta = sample_id_df
            logging.info("Using sample ID DataFrame as metadata.")
        else:
            raise ValueError("Provide metadata_file or sample_id_df.")
        sid_col = meta.columns[0]
        results = {}
        for fn in os.listdir(kraken_dir):
            if fn.endswith('_kraken_report.txt'):
                path = os.path.join(kraken_dir, fn)
                sample = fn.replace('_kraken_report.txt','')
                for line in open(path):
                    f = line.strip().split('	')
                    if len(f)<6: continue
                    perc, cov, direct, rc, taxid, sci = f[:6]
                    direct=int(direct)
                    if direct>=read_count and direct<=max_read_count and rc=='S':
                        row_meta = meta.loc[meta[sid_col]==sample].iloc[0].to_dict()
                        results[f"{sample}_{taxid}"] = {
                            'Sample': sample,'TaxID':taxid,'Scientific_Name':sci,
                            'Perc':perc,'Reads':cov,'Direct':direct,**row_meta
                        }
        out = os.path.join(kraken_dir,'merged_kraken.tsv')
        df = pd.DataFrame(results.values())
        df.to_csv(out, sep='	', index=False)
        logging.info(f"Aggregated saved to {out}")
        return out
    except Exception as e:
        logging.error(f"Error aggregating: {e}")
        return None


def generate_unfiltered_merged_tsv(kraken_dir, metadata_file=None, sample_id_df=None):
    """
    Merges all entries from Kraken summary across samples without filtering.
    """
    # Similar to aggregate but without read_count filter
    return aggregate_kraken_results(kraken_dir, metadata_file, sample_id_df, read_count=0)


def generate_abundance_plots(merged_tsv, top_N, col_filter, pat_to_keep):
    """
    Creates bar plots of abundance by category using Plotly.
    """
    try:
        df = pd.read_csv(merged_tsv, sep='	')
        df = df[df['Scientific_Name']!='Homo sapiens']
        if col_filter: df=df[~df['Scientific_Name'].isin(col_filter)]
        if pat_to_keep:df=df[df['Scientific_Name'].isin(pat_to_keep)]
        for category in df.columns:
            if category in ['Sample','TaxID','Perc','Reads','Direct']: continue
            top=df['Scientific_Name'].value_counts().nlargest(top_N).index
            plot_df=df[df['Scientific_Name'].isin(top)]
            fig=px.bar(plot_df,x='Scientific_Name',y='Direct', color=category)
            out_png=f"abundance_{category}.png"
            fig.write_image(out_png)
            logging.info(f"Plot saved {out_png}")
    except Exception as e:
        logging.error(f"Error plotting: {e}")


def process_all_ranks(kraken_dir, metadata_file=None, sample_id_df=None,
                      read_count=1, max_read_count=10**30, top_N=None, col_filter=None, pat_to_keep=None):
    """
    Runs aggregation and plotting across all taxonomic ranks.
    """
    merged = generate_unfiltered_merged_tsv(kraken_dir, metadata_file, sample_id_df)
    generate_abundance_plots(merged, top_N or 10, col_filter, pat_to_keep)
    return merged


def run_multiqc(trimmomatic_output_dir):
    try:
        subprocess.run(["multiqc", "--force", trimmomatic_output_dir], check=True)
        logging.info("MultiQC complete.")
    except Exception as e:
        logging.error(f"MultiQC error: {e}")

