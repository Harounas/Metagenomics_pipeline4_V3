#!/usr/bin/env python3
"""
kraken_abundance_pipeline.py (Version 3 - Parallel + Full Helpers)

Processes Kraken2 reports, generates abundance plots, aggregates results
(with metadata or sample IDs), and supports quality control via MultiQC.
Supports:
  - Preprocessing (fastp)
  - Optional host depletion via Bowtie2
  - Optional assembly via MetaSPAdes
  - Parallel execution of multiple samples
  - Backward-compatible helper functions (process_kraken_reports, etc.)
"""

import os
import glob
import pandas as pd
import logging
import subprocess
import threading
from concurrent.futures import ProcessPoolExecutor, as_completed
import plotly.express as px

# Local imports
from .fastp import run_fastp
from .metaspades import run_spades
from .bowtie2 import run_bowtie2
from .kraken2 import run_kraken2
from concurrent.futures import ThreadPoolExecutor, as_completed

#from multiprocessing import Semaphore


# ---------------- Logging ---------------- #
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] [PID %(process)d] %(message)s"
)

# ---------------- Core sample processing ---------------- #

def process_sample(forward, reverse, base_name, bowtie2_index, kraken_db, output_dir, threads,
                   run_bowtie=True, use_precomputed_reports=False, use_assembly=False,
                   skip_preprocessing=False, skip_existing=False,
                   assembly_semaphore=None):
    """Run the full pipeline for a single sample."""
    try:
        kraken_report = os.path.join(output_dir, f"{base_name}_kraken_report.txt")
        output_report = kraken_report

        if use_precomputed_reports:
            if not os.path.exists(kraken_report):
                raise FileNotFoundError(f"Precomputed Kraken2 report not found for {base_name}")
            return kraken_report, output_report

        # --- Preprocessing / Assembly logic --- #
        if skip_preprocessing:
            contigs_file = os.path.join(output_dir, f"{base_name}_contigs.fasta")
            if not (skip_existing and os.path.exists(contigs_file)):
                logging.info(f"[{base_name}] Running SPAdes (skip_preprocessing)")
                if assembly_semaphore:
                    with assembly_semaphore:
                        contigs_file = run_spades(forward, reverse, base_name, output_dir, threads)
                else:
                    contigs_file = run_spades(forward, reverse, base_name, output_dir, threads)
            kraken_input = contigs_file
        else:
            trimmed_forward = os.path.join(output_dir, f"{base_name}_trimmed_R1.fastq.gz")
            trimmed_reverse = os.path.join(output_dir, f"{base_name}_trimmed_R2.fastq.gz")

            if not (skip_existing and os.path.exists(trimmed_forward) and os.path.exists(trimmed_reverse)):
                logging.info(f"[{base_name}] Running fastp")
                run_fastp(forward, reverse, base_name, output_dir, threads)

            unmapped_r1, unmapped_r2 = trimmed_forward, trimmed_reverse

            if run_bowtie:
                bowtie_unmapped_r1 = os.path.join(output_dir, f"{base_name}_unmapped_1.fastq.gz")
                bowtie_unmapped_r2 = os.path.join(output_dir, f"{base_name}_unmapped_2.fastq.gz")
                if not (skip_existing and os.path.exists(bowtie_unmapped_r1) and os.path.exists(bowtie_unmapped_r2)):
                    logging.info(f"[{base_name}] Running Bowtie2")
                    run_bowtie2(trimmed_forward, trimmed_reverse, base_name, bowtie2_index, output_dir, threads)
                unmapped_r1, unmapped_r2 = bowtie_unmapped_r1, bowtie_unmapped_r2

            if use_assembly:
                contigs_file = os.path.join(output_dir, f"{base_name}_contigs.fasta")
                if not (skip_existing and os.path.exists(contigs_file)):
                    logging.info(f"[{base_name}] Running SPAdes")
                    if assembly_semaphore:
                        with assembly_semaphore:
                            contigs_file = run_spades(unmapped_r1, unmapped_r2, base_name, output_dir, threads)
                    else:
                        contigs_file = run_spades(unmapped_r1, unmapped_r2, base_name, output_dir, threads)
                kraken_input = contigs_file
            else:
                kraken_input = unmapped_r1

        # --- Kraken2 classification --- #
        if not (skip_existing and os.path.exists(kraken_report)):
            logging.info(f"[{base_name}] Running Kraken2")
            if use_assembly or skip_preprocessing:
                run_kraken2(kraken_input, None, base_name, kraken_db, output_dir, threads)
            else:
                run_kraken2(unmapped_r1, unmapped_r2, base_name, kraken_db, output_dir, threads)

        # --- Postprocessing --- #
        if os.path.exists(output_report):
            logging.info(f"[{base_name}] Splitting output report")
            process_output_report(output_report, output_dir)

        return kraken_report, output_report
    except Exception as e:
        logging.error(f"[{base_name}] Error: {e}")
        return None, None

# ---------------- Parallel wrapper ---------------- #

def find_samples(input_dir):
    """
    Find paired-end FASTQs with flexible naming conventions.
    Supports:
      - *_R1.fastq.gz / *_R2.fastq.gz
      - *_R1_001.fastq.gz / *_R2_001.fastq.gz
      - *_1.fastq.gz / *_2.fastq.gz
      - sample.1.fq.gz / sample.2.fq.gz
      - .fastq or .fq extensions
    """
    patterns = ["*.fastq.gz", "*.fq.gz", "*.fastq", "*.fq"]
    samples = []
    seen = set()

    files = []
    for pat in patterns:
        files.extend(glob.glob(os.path.join(input_dir, pat)))

    files = sorted(files)

    for f in files:
        base = os.path.basename(f)

        # Try to detect if it's R1
        if any(tag in base for tag in ["_R1", "_1.", ".1.", "_1_"]):
            # Build candidate R2 filename
            r2_candidates = []
            r2_candidates.append(f.replace("_R1", "_R2"))
            r2_candidates.append(f.replace("_1.", "_2."))
            r2_candidates.append(f.replace(".1.", ".2."))
            r2_candidates.append(f.replace("_1_", "_2_"))

            for r2 in r2_candidates:
                if os.path.exists(r2):
                    # Clean sample ID
                    sid = base
                    for suffix in [
                        "_R1.fastq.gz", "_R1.fastq", "_R1_001.fastq.gz", "_R1_001.fastq",
                        "_1.fastq.gz", "_1.fq.gz", "_1.fastq", "_1.fq",
                        ".1.fastq.gz", ".1.fq.gz", ".1.fastq", ".1.fq"
                    ]:
                        sid = sid.replace(suffix, "")
                    if sid not in seen:
                        samples.append((f, r2, sid))
                        seen.add(sid)
                    break
    return samples

def process_samples_in_parallel(samples, bowtie2_index, kraken_db, output_dir, threads,
                                run_bowtie=True, use_precomputed_reports=False, use_assembly=False,
                                skip_preprocessing=False, skip_existing=False,
                                max_workers=2, max_assemblies=1):
    """Run many samples in parallel with ProcessPoolExecutor. Limits concurrent MetaSPAdes jobs."""
    results = []
    assembly_semaphore = threading.Semaphore(max_assemblies)

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                process_sample,
                f, r, b,
                bowtie2_index, kraken_db, output_dir, threads,
                run_bowtie, use_precomputed_reports, use_assembly,
                skip_preprocessing, skip_existing, assembly_semaphore
            ): b for f, r, b in samples
        }
        for fut in as_completed(futures):
            base = futures[fut]
            try:
                res = fut.result()
                results.append(res)
                logging.info(f"[DONE] {base}")
            except Exception as e:
                logging.error(f"[FAIL] {base}: {e}")
    return results

# ---------------- Report splitting ---------------- #

def process_output_report(output_report, output_dir):
    """Splits Kraken report into domain-specific files."""
    try:
        with open(output_report, 'r') as file:
            lines = file.readlines()
        current_domain, current_rows = None, []
        for line in lines:
            cols = line.strip().split("\t")
            if len(cols) < 6: continue
            rank_code = cols[3]
            if rank_code == "D":
                if current_domain:
                    save_domain_data(current_domain, current_rows, output_dir)
                current_domain, current_rows = cols[5], [line]
            else:
                current_rows.append(line)
        if current_domain:
            save_domain_data(current_domain, current_rows, output_dir)
    except Exception as e:
        logging.error(f"Error processing output report {output_report}: {e}")

def save_domain_data(domain, rows, output_dir):
    path = os.path.join(output_dir, f"{domain.replace(' ', '')}_output_report.txt")
    with open(path, "w") as f:
        f.writelines(rows)
    logging.info(f"Saved {domain} -> {path}")

# ---------------- Aggregation & plotting ---------------- #

def aggregate_kraken_results(kraken_dir, metadata_file=None, sample_id_df=None,
                             read_count=1, max_read_count=10**30):
    """Aggregates Kraken results across samples."""
    try:
        if metadata_file:
            meta = pd.read_csv(metadata_file)
        elif sample_id_df is not None:
            meta = sample_id_df
        else:
            raise ValueError("Provide metadata_file or sample_id_df.")
        sid_col = meta.columns[0]

        results = {}
        for fn in os.listdir(kraken_dir):
            if fn.endswith("_kraken_report.txt"):
                path = os.path.join(kraken_dir, fn)
                sample = fn.replace("_kraken_report.txt", "")
                for line in open(path):
                    f = line.strip().split("\t")
                    if len(f) < 6: continue
                    perc, cov, direct, rc, taxid, sci = f[:6]
                    direct = int(direct)
                    if direct >= read_count and direct <= max_read_count and rc == "S":
                        row_meta = meta.loc[meta[sid_col] == sample].iloc[0].to_dict()
                        results[f"{sample}_{taxid}"] = {
                            "Sample": sample, "TaxID": taxid, "Scientific_Name": sci,
                            "Perc": perc, "Reads": cov, "Direct": direct, **row_meta
                        }
        out = os.path.join(kraken_dir, "merged_kraken.tsv")
        pd.DataFrame(results.values()).to_csv(out, sep="\t", index=False)
        logging.info(f"Aggregated results -> {out}")
        return out
    except Exception as e:
        logging.error(f"Error aggregating: {e}")
        return None

def generate_abundance_plots(merged_tsv, top_N=10, col_filter=None, pat_to_keep=None):
    """Generate barplots of abundant taxa."""
    try:
        df = pd.read_csv(merged_tsv, sep="\t")
        df = df[df["Scientific_Name"] != "Homo sapiens"]
        if col_filter: df = df[~df["Scientific_Name"].isin(col_filter)]
        if pat_to_keep: df = df[df["Scientific_Name"].isin(pat_to_keep)]
        top = df["Scientific_Name"].value_counts().nlargest(top_N).index
        plot_df = df[df["Scientific_Name"].isin(top)]
        fig = px.bar(plot_df, x="Scientific_Name", y="Direct", color="Sample")
        fig.write_image("abundance.png")
        logging.info("Saved abundance.png")
    except Exception as e:
        logging.error(f"Plotting error: {e}")

# ---------------- MultiQC ---------------- #

def run_multiqc(out_dir):
    try:
        subprocess.run(["multiqc", "--force", out_dir], check=True)
        logging.info("MultiQC complete.")
    except Exception as e:
        logging.error(f"MultiQC error: {e}")

# ---------------- Backward compatibility helpers ---------------- #

def process_kraken_reports(kraken_dir):
    """Process all Kraken reports in a directory into domain-specific files."""
    for fn in os.listdir(kraken_dir):
        if fn.endswith("_kraken_report.txt"):
            process_output_report(os.path.join(kraken_dir, fn), kraken_dir)

def process_output_reports(kraken_dir):
    """Process all *_output_report.txt files in a directory."""
    for fn in os.listdir(kraken_dir):
        if fn.endswith("_output_report.txt"):
            process_output_report(os.path.join(kraken_dir, fn), kraken_dir)

def generate_unfiltered_merged_tsv(kraken_dir, metadata_file=None, sample_id_df=None):
    """Aggregate Kraken results across all samples with no read_count filtering."""
    return aggregate_kraken_results(kraken_dir, metadata_file, sample_id_df, read_count=0)

def process_all_ranks(kraken_dir, metadata_file=None, sample_id_df=None,
                      read_count=1, max_read_count=10**30, top_N=None,
                      col_filter=None, pat_to_keep=None):
    """Aggregate results across all ranks and generate abundance plots."""
    merged = generate_unfiltered_merged_tsv(kraken_dir, metadata_file, sample_id_df)
    generate_abundance_plots(merged, top_N or 10, col_filter, pat_to_keep)
    return merged
