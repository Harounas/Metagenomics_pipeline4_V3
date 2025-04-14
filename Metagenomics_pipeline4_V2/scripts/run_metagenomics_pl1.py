#!/usr/bin/env python3

import os
import glob
import argparse
import pandas as pd
import sys
import logging
from Metagenomics_pipeline4_V2.kraken_abundance_pipeline import (
    process_sample, 
    aggregate_kraken_results, 
    generate_abundance_plots,
    run_multiqc,
    process_kraken_reports
)

# Configure logging
logging.basicConfig(
    level=logging.INFO, 
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("metagenomics_pipeline.log"),
        logging.StreamHandler()
    ]
)

def create_sample_id_df(input_dir):
    sample_ids = []
    for f in glob.glob(os.path.join(input_dir, "*_R1*.fastq*")):
        sample_id = os.path.basename(f)
        for pattern in ["_R1_001.fastq.gz", "_R1_001.fastq", "_R1.fastq.gz", "_R1.fastq", "R1.fastq.gz", "R1.fastq", "_R1_001", "_R1"]:
            sample_id = sample_id.replace(pattern, "")
        sample_ids.append(sample_id)
    return pd.DataFrame(sample_ids, columns=["Sample_IDs"])


def validate_inputs(args):
    if not os.path.isdir(args.input_dir):
        logging.error(f"Input directory '{args.input_dir}' not found.")
        sys.exit(1)
    if not os.path.isdir(args.kraken_db):
        logging.error(f"Kraken database directory '{args.kraken_db}' not found.")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Metagenomics pipeline")
    parser.add_argument("--input_dir", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--kraken_db", required=True)
    parser.add_argument("--bowtie2_index", required=False)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--no_metadata", action='store_true')
    parser.add_argument("--use_assembly", action='store_true')
    parser.add_argument("--paired_assembly", action='store_true')
    parser.add_argument("--skip_preprocessing", action='store_true')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    validate_inputs(args)

    sample_id_df = create_sample_id_df(args.input_dir) if args.no_metadata else None

    for forward in glob.glob(os.path.join(args.input_dir, "*_R1*.fastq*")):
        base_name = os.path.basename(forward)
        for pattern in ["_R1_001.fastq.gz", "_R1.fastq.gz", "_R1.fastq", "_R1"]:
            base_name = base_name.replace(pattern, "")

        reverse = None
        if args.paired_assembly or not args.use_assembly:
            reverse_candidates = [
                os.path.join(args.input_dir, f"{base_name}_R2_001.fastq.gz"),
                os.path.join(args.input_dir, f"{base_name}_R2.fastq.gz"),
                os.path.join(args.input_dir, f"{base_name}_R2.fastq")
            ]
            reverse = next((f for f in reverse_candidates if os.path.isfile(f)), None)

        process_sample(
            forward=forward,
            reverse=reverse,
            base_name=base_name,
            bowtie2_index=args.bowtie2_index,
            kraken_db=args.kraken_db,
            output_dir=args.output_dir,
            threads=args.threads,
            run_bowtie=bool(args.bowtie2_index),
            use_precomputed_reports=False,
            use_assembly=args.use_assembly,
            skip_preprocessing=args.skip_preprocessing,
            skip_existing=False
        )

    merged_tsv = aggregate_kraken_results(args.output_dir, sample_id_df=sample_id_df)

    process_kraken_reports(args.output_dir)

    generate_abundance_plots(merged_tsv)

    run_multiqc(args.output_dir)


if __name__ == "__main__":
    main()
