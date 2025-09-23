#!/usr/bin/env python3
"""
run_metagenomics_pl2.py

Metagenomics pipeline including Kraken2 classification, optional assembly-based contig
extraction, geNomad detection, clustering, and Diamond annotation with downstream analyses.
Supports:
  - Preprocessing (fastp)
  - Optional host depletion via Bowtie2
  - Optional assembly via MetaSPAdes
  - Parallel execution of multiple samples with controlled SPAdes concurrency
"""

import os
import glob
import argparse
import pandas as pd
import sys
import logging
import subprocess
import csv
from Bio import SeqIO, Entrez

# Import functions from pipeline modules
from Metagenomics_pipeline4_V2.kraken_abundance_pipeline import (
    process_samples_in_parallel,
    find_samples,
    aggregate_kraken_results,
    generate_abundance_plots,
    run_multiqc,
    process_kraken_reports,
    generate_unfiltered_merged_tsv,
    process_all_ranks,
    process_output_reports
)
from Metagenomics_pipeline4_V2.ref_based_assembly import ref_based
from Metagenomics_pipeline4_V2.deno_ref_assembly2 import deno_ref_based
from Metagenomics_pipeline4_V2 import extract_contigs_diamond
from Metagenomics_pipeline4_V2.alignment_summary import run_alignment_summary
from Metagenomics_pipeline4_V2.extract_contigs_diamond import process_virus_contigs
from Metagenomics_pipeline4_V2.process_clustered_contigs import process_clustered_contigs

# ---------------- Logging ---------------- #
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("metagenomics_pipeline.log"),
        logging.StreamHandler()
    ]
)

# ---------------- Sample utilities ---------------- #

def create_sample_id_df(input_dir):
    """Create DataFrame of sample IDs based on FASTQ names."""
    sample_ids = []
    for f in glob.glob(os.path.join(input_dir, "*_R1*.fastq*")):
        sid = os.path.basename(f)
        for pat in [
            "_R1_001.fastq.gz", "_R1_001.fastq", "_R1.fastq.gz",
            "_R1.fastq", "R1.fastq.gz", "R1.fastq", "_R1_001", "_R1"
        ]:
            sid = sid.replace(pat, "")
        sample_ids.append(sid)
    return pd.DataFrame(sample_ids, columns=["Sample_IDs"])

def validate_inputs(args):
    """Validate paths and required arguments."""
    if not os.path.isdir(args.input_dir):
        logging.error(f"Input directory '{args.input_dir}' not found.")
        sys.exit(1)
    if not os.path.isdir(args.kraken_db):
        logging.error(f"Kraken DB '{args.kraken_db}' not found.")
        sys.exit(1)
    if args.bowtie2_index and not os.path.exists(args.bowtie2_index + ".1.bt2"):
        logging.error(f"Bowtie2 index '{args.bowtie2_index}' not found.")
        sys.exit(1)
    if args.metadata_file and not os.path.isfile(args.metadata_file):
        logging.error(f"Metadata file '{args.metadata_file}' not found.")
        sys.exit(1)
    if args.use_precomputed_reports and not glob.glob(os.path.join(args.output_dir, "*_report.txt")):
        logging.error("No precomputed Kraken reports found in output directory")
        sys.exit(1)
    if args.run_genomad and not args.genomad_db:
        logging.error("You must provide --genomad_db when using --run_genomad")
        sys.exit(1)

# ---------------- Parallel processing ---------------- #

def process_samples(args):
    """Process all samples in parallel with controlled SPAdes concurrency."""
    run_bowtie = not args.no_bowtie2 and args.bowtie2_index is not None

    samples = find_samples(args.input_dir)
    if not samples:
        logging.error("❌ No samples found in input directory.")
        sys.exit(1)

    logging.info(f"Found {len(samples)} samples. Running {args.parallel} jobs in parallel "
                 f"(max {args.max_assemblies} SPAdes assemblies at once).")

    process_samples_in_parallel(
        samples,
        bowtie2_index=args.bowtie2_index,
        kraken_db=args.kraken_db,
        output_dir=args.output_dir,
        threads=args.threads,
        run_bowtie=run_bowtie,
        use_precomputed_reports=args.use_precomputed_reports,
        use_assembly=args.use_assembly,
        skip_preprocessing=args.skip_preprocessing,
        skip_existing=args.skip_existing,
        max_workers=args.parallel,
        max_assemblies=args.max_assemblies
    )

# ---------------- Metadata ---------------- #

def handle_metadata(args):
    """Aggregate Kraken reports and merge with metadata or sample IDs."""
    if args.no_metadata:
        df = create_sample_id_df(args.input_dir)
        df.to_csv(os.path.join(args.output_dir, "sample_ids.csv"), index=False)
        return aggregate_kraken_results(
            args.output_dir,
            sample_id_df=df,
            read_count=args.read_count,
            max_read_count=args.max_read_count
        )
    return aggregate_kraken_results(
        args.output_dir,
        metadata_file=args.metadata_file,
        read_count=args.read_count,
        max_read_count=args.max_read_count
    )

# ---------------- Main ---------------- #

def main():
    parser = argparse.ArgumentParser(
        description="Metagenomics pipeline for taxonomic classification and analysis"
    )
    # Core
    parser.add_argument("--kraken_db", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--input_dir", required=True)
    parser.add_argument("--bowtie2_index")
    parser.add_argument("--threads", type=int, default=8,
                        help="Threads per sample job (fastp/Bowtie2/SPAdes/Kraken2)")
    parser.add_argument("--parallel", type=int, default=2,
                        help="Number of samples processed in parallel")
    parser.add_argument("--max_assemblies", type=int, default=1,
                        help="Max concurrent MetaSPAdes assemblies")
    # Kraken/metadata
    parser.add_argument("--metadata_file")
    parser.add_argument("--read_count", type=int, default=1)
    parser.add_argument("--max_read_count", type=int, default=5000000000)
    parser.add_argument("--top_N", type=int, default=10)
    parser.add_argument("--no_bowtie2", action="store_true")
    parser.add_argument("--no_metadata", action="store_true")
    parser.add_argument("--use_precomputed_reports", action="store_true")
    parser.add_argument("--use_assembly", action="store_true")
    parser.add_argument("--skip_preprocessing", action="store_true")
    parser.add_argument("--skip_existing", action="store_true")
    parser.add_argument("--skip_multiqc", action="store_true")
    parser.add_argument("--skip_reports", action="store_true")
    parser.add_argument("--process_all_ranks", action="store_true")
    parser.add_argument("--col_filter", type=str, nargs="+")
    parser.add_argument("--pat_to_keep", type=str, nargs="+")
    # geNomad / Diamond
    parser.add_argument("--diamond", action="store_true")
    parser.add_argument("--diamond_db")
    parser.add_argument("--nr_path", type=str, help="Path to nr FASTA file containing virus accession and name")
    parser.add_argument("--skip_diamond", action="store_true")
    parser.add_argument("--run_genomad", action="store_true")
    parser.add_argument("--genomad_db")
    parser.add_argument("--skip_genomad", action="store_true")
    # Downstream
    parser.add_argument("--run_alignment", action="store_true")
    parser.add_argument("--max_workers", type=int, default=5)
    parser.add_argument("--bwa_threads", type=int, default=4)
    parser.add_argument("--run_scaffolding", action="store_true")
    parser.add_argument("--run_ref_base", action="store_true")
    parser.add_argument("--run_deno_ref", action="store_true")

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    validate_inputs(args)

    # 🔹 Parallel sample execution
    process_samples(args)
    merged_tsv = handle_metadata(args)

    if not args.skip_reports:
        logging.info("Processing Kraken reports…")
        process_kraken_reports(args.output_dir)
        logging.info("Processing output reports…")
        process_output_reports(args.output_dir)

    # ---------------- Diamond + geNomad ---------------- #
    if args.diamond and not args.skip_diamond:
        if args.use_assembly:
            logging.info("Running contig extraction and Diamond annotation…")
            long_contigs_fasta = extract_contigs_diamond.extract_long_contigs_kraken(
                base_contigs_dir=args.output_dir,
                output_tsv=os.path.join(args.output_dir, "long_contigs_summary.tsv")
            )
            if args.run_genomad and not args.skip_genomad:
                genomad_input_fasta = os.path.join(args.output_dir, "merged_contigs_genomad.fasta")
                genomad_out_dir = os.path.join(args.output_dir, "genomad_output")
                clustered_out_dir = os.path.join(args.output_dir, "clustered_output")
                final_long_clustered_fasta = os.path.join(args.output_dir, "clustered_long_contigs.fasta")

                extract_contigs_diamond.extract_and_merge_contigs_genomad(
                    base_contigs_dir=args.output_dir,
                    output_fasta=genomad_input_fasta
                )
                extract_contigs_diamond.run_genomad(
                    input_fasta=genomad_input_fasta,
                    output_dir=genomad_out_dir,
                    genomad_db=args.genomad_db,
                    threads=args.threads
                )
                clustered_fasta = extract_contigs_diamond.cluster_contigs(
                    virus_fasta=genomad_input_fasta,
                    output_dir=clustered_out_dir,
                    threads=args.threads
                )
                extract_contigs_diamond.extract_long_contigs(
                    input_fasta=clustered_fasta,
                    output_fasta=final_long_clustered_fasta
                )
                extract_contigs_diamond.run_diamond(
                    diamond_db=args.diamond_db,
                    query_file=final_long_clustered_fasta,
                    output_file=os.path.join(args.output_dir, "results_clustered.m8"),
                    threads=args.threads
                )
                process_virus_contigs(
                    fasta_file=args.nr_path,
                    diamond_results_file=os.path.join(args.output_dir, "results_clustered.m8"),
                    output_dir=args.output_dir
                )
        else:
            logging.warning("Diamond requested but --use_assembly not provided. Skipping contig analysis.")

    if not args.skip_multiqc:
        run_multiqc(args.output_dir)

if __name__ == "__main__":
    main()

