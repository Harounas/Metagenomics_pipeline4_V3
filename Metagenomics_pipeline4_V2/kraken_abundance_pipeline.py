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

        if use_precomputed_reports:
            if not os.path.exists(kraken_report):
                raise FileNotFoundError(f"Precomputed Kraken2 report not found: {kraken_report}")
            return kraken_report

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

        return kraken_report

    except Exception as e:
        logging.error(f"Error processing sample {base_name}: {e}")
        return None


def run_multiqc(output_dir):
    try:
        subprocess.run(["multiqc", output_dir], check=True)
        logging.info("MultiQC report generated successfully.")
    except Exception as e:
        logging.error(f"Error running MultiQC: {e}")
