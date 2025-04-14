import os
import glob
import argparse
import pandas as pd
import sys
import logging
import subprocess
from Metagenomics_pipeline4_V2.kraken_abundance_pipeline import (
    process_sample, 
    aggregate_kraken_results, 
    generate_abundance_plots,
    run_multiqc,
    process_kraken_reports,
    generate_unfiltered_merged_tsv,
    process_all_ranks
)
from Metagenomics_pipeline4_V2.ref_based_assembly import ref_based
from Metagenomics_pipeline4_V2.deno_ref_assembly2 import deno_ref_based

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
    required_paths = [args.input_dir, args.kraken_db]
    if args.bowtie2_index:
        required_paths.append(args.bowtie2_index + ".1.bt2")
    if args.metadata_file:
        required_paths.append(args.metadata_file)

    for path in required_paths:
        if not os.path.exists(path):
            logging.error(f"Required path '{path}' not found.")
            sys.exit(1)

    if args.use_precomputed_reports and not glob.glob(os.path.join(args.output_dir, "*_report.txt")):
        logging.error("No precomputed Kraken reports found in output directory")
        sys.exit(1)

def process_samples(args):
    run_bowtie = not args.no_bowtie2 and args.bowtie2_index is not None
    for forward in glob.glob(os.path.join(args.input_dir, "*_R1*.fastq*")):
        base_name = os.path.basename(forward)
        for pattern in ["_R1_001.fastq.gz", "_R1_001.fastq", "_R1.fastq.gz", "_R1.fastq", "R1.fastq.gz", "R1.fastq", "_R1_001", "_R1"]:
            base_name = base_name.replace(pattern, "")

        reverse = None
        if not args.use_assembly or args.paired_assembly:
            reverse_candidates = [
                os.path.join(args.input_dir, f"{base_name}_R2_001.fastq.gz"),
                os.path.join(args.input_dir, f"{base_name}_R2.fastq.gz"),
                os.path.join(args.input_dir, f"{base_name}_R2.fastq")
            ]
            reverse = next((f for f in reverse_candidates if os.path.isfile(f)), None)
            if not reverse and not args.use_assembly:
                logging.warning(f"No matching R2 file found for {base_name}. Skipping.")
                continue

        logging.info(f"Processing sample {base_name} (Assembly: {args.use_assembly})")
        process_sample(
            forward=forward,
            reverse=reverse,
            base_name=base_name,
            bowtie2_index=args.bowtie2_index,
            kraken_db=args.kraken_db,
            output_dir=args.output_dir,
            threads=args.threads,
            run_bowtie=run_bowtie,
            use_precomputed_reports=args.use_precomputed_reports,
            use_assembly=args.use_assembly,
            skip_preprocessing=args.skip_preprocessing,
            skip_existing=args.skip_existing
        )

def handle_metadata(args):
    if args.no_metadata:
        sample_id_df = create_sample_id_df(args.input_dir)
        sample_id_df.to_csv(os.path.join(args.output_dir, "sample_ids.csv"), index=False)
        return aggregate_kraken_results(args.output_dir, sample_id_df=sample_id_df, read_count=args.read_count, max_read_count=args.max_read_count)
    return aggregate_kraken_results(args.output_dir, metadata_file=args.metadata_file, read_count=args.read_count, max_read_count=args.max_read_count)

def main():
    parser = argparse.ArgumentParser(description="Metagenomics pipeline for taxonomic classification and analysis", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # arguments definition here (unchanged)...

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    validate_inputs(args)

    process_samples(args)
    merged_tsv_path = handle_metadata(args)

    if not args.skip_multiqc:
        run_multiqc(args.output_dir)

    if not args.skip_reports:
        process_kraken_reports(args.output_dir)
        domains = ["Bacteria", "Viruses", "Archaea", "Eukaryota"]
        domain_flags = [args.bacteria, args.virus, args.archaea, args.eukaryota]

        for domain, flag in zip(domains, domain_flags):
            if flag:
                merged_tsv = aggregate_kraken_results(args.output_dir, args.metadata_file, None, args.read_count, args.max_read_count)
                if merged_tsv:
                    generate_abundance_plots(merged_tsv, args.top_N, args.col_filter, args.pat_to_keep)

if __name__ == "__main__":
    main()
