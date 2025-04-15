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


def process_filename(fname, token="kraken", new_ext=".txt"):
    base, _ = os.path.splitext(fname)
    return "_".join(t for t in base.split("_") if t.lower() != token.lower()) + new_ext

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
                    domain_output_filename=process_filename(domain_output_filename)
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

def aggregate_kraken_results(kraken_dir, metadata_file=None, sample_id_df=None,
                             read_count=1, max_read_count=10**30):
    """
    Aggregates Kraken results (species-level) and merges metadata (or sample IDs if metadata_file is None)
    into a single TSV file.

    Parameters:
      kraken_dir (str): Directory containing Kraken report files.
      metadata_file (str, optional): Path to metadata CSV.
      sample_id_df (DataFrame, optional): DataFrame of sample IDs.
      read_count (int): Minimum read count threshold.
      max_read_count (int): Maximum read count threshold.

    Returns:
      str: Path to the generated merged TSV file.
    """
    try:
        # Load metadata
        if metadata_file:
            metadata = pd.read_csv(metadata_file, sep=",")
            logging.info("Using metadata from the provided file.")
        elif sample_id_df is not None:
            metadata = sample_id_df
            logging.info("Using sample IDs as metadata.")
        else:
            raise ValueError("Either metadata_file or sample_id_df must be provided.")

        sample_id_col = metadata.columns[0]
        aggregated_results = {}

        for file_name in os.listdir(kraken_dir):
            if file_name.endswith("_report.txt"):
                with open(os.path.join(kraken_dir, file_name), 'r') as f:
                    for line in f:
                        fields = line.strip().split('\t')
                        if len(fields) < 6:
                            continue
                        perc_frag_cover = fields[0]
                        nr_frag_cover = fields[1]
                        nr_frag_direct_at_taxon = int(fields[2])
                        rank_code_field = fields[3]
                        ncbi_ID = fields[4]
                        scientific_name = fields[5]
                        parts = file_name.split('_')
                        extracted_part = '_'.join(parts[:-2])
                        sampleandtaxonid = f"{extracted_part}{ncbi_ID}"

                        if rank_code_field in ['S', 'S1', 'S2', 'S3'] and (read_count <= nr_frag_direct_at_taxon <= max_read_count):
                            if extracted_part in metadata[sample_id_col].unique():
                                sample_metadata = metadata.loc[metadata[sample_id_col] == extracted_part].iloc[0].to_dict()
                                aggregated_results[sampleandtaxonid] = {
                                    'Perc_frag_cover': perc_frag_cover,
                                    'Nr_frag_cover': nr_frag_cover,
                                    'Nr_frag_direct_at_taxon': nr_frag_direct_at_taxon,
                                    'Rank_code': rank_code_field,
                                    'NCBI_ID': ncbi_ID,
                                    'Scientific_name': scientific_name,
                                    'SampleID': extracted_part,
                                    **sample_metadata
                                }

        merged_tsv_path = os.path.join(kraken_dir, "merged_kraken.tsv")
        with open(merged_tsv_path, 'w') as f:
            headers = ['Perc_frag_cover', 'Nr_frag_cover', 'Nr_frag_direct_at_taxon',
                       'Rank_code', 'NCBI_ID', 'Scientific_name', 'SampleID'] + metadata.columns[1:].tolist()
            f.write("\t".join(headers) + "\n")
            for data in aggregated_results.values():
                f.write("\t".join(str(data[col]) for col in headers) + "\n")

        logging.info(f"Aggregated results saved to {merged_tsv_path}")
        return merged_tsv_path

    except Exception as e:
        logging.error(f"Error aggregating Kraken results: {e}")
        return None


def generate_unfiltered_merged_tsv(kraken_dir, metadata_file=None, sample_id_df=None):
    """
    Generates an unfiltered merged TSV file containing all Kraken report data from the specified directory.

    Parameters:
      kraken_dir (str): Directory with Kraken report files.
      metadata_file (str, optional): Path to metadata CSV.
      sample_id_df (DataFrame, optional): DataFrame of sample IDs if metadata_file is not provided.

    Returns:
      str: Path to the unfiltered merged TSV file.
    """
    try:
        if metadata_file:
            metadata = pd.read_csv(metadata_file, sep=",")
            logging.info("Using metadata from the provided file.")
        elif sample_id_df is not None:
            metadata = sample_id_df
            logging.info("Using sample IDs as metadata.")
        else:
            raise ValueError("Either metadata_file or sample_id_df must be provided.")

        sample_id_col = metadata.columns[0]
        unfiltered_results = {}

        for file_name in os.listdir(kraken_dir):
            if file_name.endswith("_report.txt"):
                with open(os.path.join(kraken_dir, file_name), 'r') as f:
                    for line in f:
                        fields = line.strip().split('\t')
                        if len(fields) < 6:
                            continue
                        perc_frag_cover = fields[0]
                        nr_frag_cover = fields[1]
                        nr_frag_direct_at_taxon = int(fields[2])
                        rank_code_field = fields[3]
                        ncbi_ID = fields[4]
                        scientific_name = fields[5]
                        parts = file_name.split('_')
                        extracted_part = '_'.join(parts[:-2])
                        sampleandtaxonid = f"{extracted_part}{ncbi_ID}"

                        if extracted_part in metadata[sample_id_col].unique():
                            sample_metadata = metadata.loc[metadata[sample_id_col] == extracted_part].iloc[0].to_dict()
                            unfiltered_results[sampleandtaxonid] = {
                                'Perc_frag_cover': perc_frag_cover,
                                'Nr_frag_cover': nr_frag_cover,
                                'Nr_frag_direct_at_taxon': nr_frag_direct_at_taxon,
                                'Rank_code': rank_code_field,
                                'NCBI_ID': ncbi_ID,
                                'Scientific_name': scientific_name,
                                'SampleID': extracted_part,
                                **sample_metadata
                            }

        merged_tsv_path = os.path.join(kraken_dir, "merged_kraken_all_ranks_unfiltered.tsv")
        with open(merged_tsv_path, 'w') as f:
            headers = ['Perc_frag_cover', 'Nr_frag_cover', 'Nr_frag_direct_at_taxon',
                       'Rank_code', 'NCBI_ID', 'Scientific_name', 'SampleID'] + metadata.columns[1:].tolist()
            f.write("\t".join(headers) + "\n")
            for data in unfiltered_results.values():
                f.write("\t".join(str(data[col]) for col in headers) + "\n")

        logging.info(f"Unfiltered merged Kraken results saved to {merged_tsv_path}")
        return merged_tsv_path

    except Exception as e:
        logging.error(f"Error generating unfiltered merged TSV: {e}")
        return None


def generate_abundance_plots(merged_tsv_path, top_N, col_filter, pat_to_keep):
    """
    Generates abundance plots for viral and bacterial data from an aggregated Kraken TSV.
    Applies optional filtering (removing human reads, and user-specified taxa).
    
    Parameters:
      merged_tsv_path (str): Path to the merged Kraken TSV.
      top_N (int): Limit to top N categories.
      col_filter (list): List of taxa names to remove.
      pat_to_keep (list): List of taxa names to exclusively retain.
    """
    try:
        df = pd.read_csv(merged_tsv_path, sep="\t")
        # Clean column names and strip string values
        df.columns = df.columns.str.replace('/', '_').str.replace(' ', '_')
        df = df.apply(lambda col: col.map(lambda x: x.strip() if isinstance(x, str) else x))
        # Remove human reads
        df = df[df['Scientific_name'] != 'Homo sapiens']
        if col_filter:
            df = df[~df['Scientific_name'].isin(col_filter)]
        if pat_to_keep:
            df = df[df['Scientific_name'].isin(pat_to_keep)]

        # Generate plots for both viral and bacterial abundance
        for focus, filter_str, plot_title in [
            ('Virus_Type', 'Virus', 'Viral'),
            ('Bacteria_Type', 'Virus', 'Bacterial')
        ]:
            if focus == 'Bacteria_Type':
                df_focus = df[~df['Scientific_name'].str.contains(filter_str, case=False, na=False)]
            else:
                df_focus = df[df['Scientific_name'].str.contains(filter_str, case=False, na=False)]
            df_focus = df_focus.rename(columns={'Scientific_name': focus})

            if top_N:
                top_N_categories = df_focus[focus].value_counts().head(top_N).index
                df_focus = df_focus[df_focus[focus].isin(top_N_categories)]

            categorical_cols = df_focus.select_dtypes(include=['object']).columns.tolist()
            if focus in categorical_cols:
                categorical_cols.remove(focus)

            for col in categorical_cols:
                grouped_sum = df_focus.groupby([focus, col])['Nr_frag_direct_at_taxon'].mean().reset_index()
                # Define color mapping: use preset colors if possible, else generate random ones.
                colordict = defaultdict(int)
                preset_colors = ['#000000','#FF0000','#556B2F','#ADD8E6','#6495ED','#00FF00',
                                 '#0000FF','#FFFF00','#00FFFF','#FF00FF','#C0C0C0','#808080',
                                 '#800000','#808000','#008000','#008080','#000080','#CD5C5C',
                                 '#DAA520','#FFA500','#F0E68C','#ADFF2F','#2F4F4F','#E0FFFF',
                                 '#4169E1','#8A2BE2','#4B0082','#EE82EE','#D2691E','#BC8F8F',
                                 '#800080','#DDA0DD','#FF1493','#8B4513','#A0522D','#708090',
                                 '#B0C4DE','#FFFFF0','#DCDCDC','#FFEFD5','#F5DEB3','#7FFFD4',
                                 '#FFC0CB','#A52A2A','#040720','#34282C','#3B3131','#3A3B3C',
                                 '#52595D','#FFFFFF','#FFFFF4','#FFF9E3']
                unique_targets = grouped_sum[focus].unique()
                if len(unique_targets) <= len(preset_colors):
                    for target, color in zip(unique_targets, preset_colors[:len(unique_targets)]):
                        colordict[target] = color
                else:
                    random_colors = [f"#{random.randint(0, 0xFFFFFF):06X}" for _ in range(len(unique_targets))]
                    for target, color in zip(unique_targets, random_colors):
                        colordict[target] = color

                # Dynamic plot dimensions and font size
                plot_width = 1100 + 5 * len(grouped_sum[col].unique())
                plot_height = 800 + 5 * len(grouped_sum[col].unique())
                font_size = max(10, 14 - len(grouped_sum[col].unique()) // 10)

                fig = px.bar(
                    grouped_sum,
                    x=col,
                    y='Nr_frag_direct_at_taxon',
                    color=focus,
                    color_discrete_map=colordict,
                    title=f"{plot_title} Abundance by {col}"
                )
                summary_csv_path = os.path.join(f"{plot_title}_summary.csv")
                grouped_sum.to_csv(summary_csv_path, index=False)
                fig.update_layout(
                    xaxis=dict(tickfont=dict(size=font_size), tickangle=45),
                    yaxis=dict(tickfont=dict(size=font_size)),
                    title=dict(text=f'Average {plot_title} Abundance by {col}', x=0.5, font=dict(size=16)),
                    bargap=0.5,
                    legend=dict(
                        font=dict(size=font_size),
                        x=1, y=1,
                        traceorder='normal',
                        orientation='v',
                        itemwidth=30,
                        itemsizing='constant',
                        itemclick='toggleothers',
                        itemdoubleclick='toggle'
                    ),
                    width=plot_width,
                    height=plot_height
                )
                out_img = f"{plot_title}_Abundance_by_{col}.png"
                fig.write_image(out_img, format='png', scale=3, width=1920, height=1080)
                logging.info(f"Abundance plot saved to {out_img}")

    except Exception as e:
        logging.error(f"Error generating abundance plots: {e}")

def process_all_ranks(kraken_dir, metadata_file=None, sample_id_df=None,
                      read_count=1, max_read_count=10**30, top_N=None, col_filter=None, pat_to_keep=None):
    """
    Processes Kraken results by generating abundance plots for multiple rank codes (S, K, G, F, D)
    and creates an unfiltered merged TSV.

    Parameters:
      kraken_dir (str): Directory with Kraken report files.
      metadata_file (str, optional): Metadata CSV file path.
      sample_id_df (DataFrame, optional): DataFrame of sample IDs (if metadata not provided).
      read_count (int): Minimum read count.
      max_read_count (int): Maximum read count.
      top_N (int): Top N categories to plot.
      col_filter (list): List of taxa to filter out.
      pat_to_keep (list): List of taxa to exclusively keep.

    Returns:
      str: Path to the unfiltered merged TSV.
    """
    unfiltered_tsv = generate_unfiltered_merged_tsv(kraken_dir, metadata_file, sample_id_df)
    rank_codes = ['S', 'K', 'G', 'F', 'D']
    domain_rank_mapping = {
        "Bacteria": ["D", "F", "K", "G", "S"],
        "Eukaryota": ["D", "F", "K", "G", "S"],
        "Viruses": ["D", "F", "K", "G", "S"]
    }

    for domain, ranks in domain_rank_mapping.items():
        for rank in ranks:
            # Generate merged Kraken results for each rank and domain
            merged_tsv = aggregate_kraken_results(kraken_dir, metadata_file, sample_id_df,
                                                  read_count, max_read_count)
            if merged_tsv:
                # Generate abundance plots for each rank and domain
                generate_abundance_plots(merged_tsv, top_N, col_filter, pat_to_keep)
                
                # Save the aggregated results to specific domain rank files
                output_filename = f"merged_kraken_{rank}_{domain}.tsv"
                output_path = os.path.join(kraken_dir, output_filename)
                logging.info(f"Saving merged Kraken results for {domain} at rank {rank} to {output_path}")
                # Save the merged results for each domain/rank (add your save logic here)
                # Example: merged_tsv.to_csv(output_path, sep="\t")

    return unfiltered_tsv

def run_multiqc(trimmomatic_output_dir):
    """
    Runs MultiQC on all files in the specified directory.

    Parameters:
      trimmomatic_output_dir (str): Directory containing files to summarize with MultiQC.
    """
    try:
        subprocess.run(["multiqc", trimmomatic_output_dir], check=True)
        logging.info("MultiQC report generated successfully.")
    except Exception as e:
        logging.error(f"Error running MultiQC: {e}")
