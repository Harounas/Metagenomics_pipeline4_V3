#!/usr/bin/env python3
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
    """
    Create a DataFrame with sample IDs based on the input FASTQ file names.
    """
    sample_ids = []
    for f in glob.glob(os.path.join(input_dir, "*_R1*.fastq*")):
        sample_id = os.path.basename(f)
        for pattern in [
            "_R1_001.fastq.gz", "_R1_001.fastq", 
            "_R1.fastq.gz", "_R1.fastq", 
            "R1.fastq.gz", "R1.fastq",
            "_R1_001", "_R1"
        ]:
            sample_id = sample_id.replace(pattern, "")
        sample_ids.append(sample_id)
    sample_id_df = pd.DataFrame(sample_ids, columns=["Sample_IDs"])
    return sample_id_df

def validate_inputs(args):
    """Validate input arguments and existence of required files and directories."""
    if not os.path.isdir(args.input_dir):
        logging.error(f"Input directory '{args.input_dir}' not found.")
        sys.exit(1)
    if not os.path.isdir(args.kraken_db):
        logging.error(f"Kraken database directory '{args.kraken_db}' not found.")
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

def process_samples(args):
    """Process all samples according to the specified parameters."""
    run_bowtie = not args.no_bowtie2 and args.bowtie2_index is not None

    # Loop over all forward FASTQ files matching a pattern.
    for forward in glob.glob(os.path.join(args.input_dir, "*_R1*.fastq*")):
        base_name = os.path.basename(forward)
        for pattern in [
            "_R1_001.fastq.gz", "_R1_001.fastq", 
            "_R1.fastq.gz", "_R1.fastq", 
            "R1.fastq.gz", "R1.fastq",
            "_R1_001", "_R1"
        ]:
            base_name = base_name.replace(pattern, "")
        
        # Determine matching reverse file (if available)
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
            use_precomputed_reports=args.use_precomputed_reports
        )

def handle_metadata(args):
    """
    Process metadata input. If no metadata file is provided (--no_metadata flag),
    create a CSV using sample IDs extracted from input FASTQ filenames and use that.
    Then, aggregate Kraken results using the metadata.
    """
    if args.no_metadata:
        sample_id_df = create_sample_id_df(args.input_dir)
        logging.info("Using sample IDs as metadata.")
        sample_id_df.to_csv(os.path.join(args.output_dir, "sample_ids.csv"), index=False)
        return aggregate_kraken_results(
            args.output_dir, 
            sample_id_df=sample_id_df, 
            read_count=args.read_count, 
            max_read_count=args.max_read_count
        )
    else:
        return aggregate_kraken_results(
            args.output_dir, 
            metadata_file=args.metadata_file, 
            read_count=args.read_count, 
            max_read_count=args.max_read_count
        )

def main():
    parser = argparse.ArgumentParser(
        description="Metagenomics pipeline for taxonomic classification and analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # Required arguments
    parser.add_argument("--kraken_db", required=True, help="Path to Kraken2 database.")
    parser.add_argument("--output_dir", required=True, help="Directory to save output files.")
    parser.add_argument("--input_dir", required=True, help="Directory containing input FASTQ files.")
    # Optional arguments
    parser.add_argument("--bowtie2_index", help="Path to Bowtie2 index (for host depletion).")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads to use.")
    parser.add_argument("--metadata_file", help="Path to the metadata CSV file.")
    parser.add_argument("--read_count", type=int, default=1, help="Minimum read count threshold.")
    parser.add_argument("--top_N", type=int, default=10000, help="Show top N most abundant taxa.")
    parser.add_argument("--max_read_count", type=int, default=5000000000, help="Maximum read count threshold.")
    # Flags for primary processing
    parser.add_argument("--no_bowtie2", action='store_true', help="Skip Bowtie2 host depletion.")
    parser.add_argument("--no_metadata", action='store_true', help="Use sample IDs instead of metadata file.")
    parser.add_argument("--use_precomputed_reports", action='store_true', help="Use existing Kraken reports instead of running Kraken2.")
    parser.add_argument("--use_assembly", action='store_true', help="Perform de novo assembly before classification.")
    parser.add_argument("--paired_assembly", action='store_true', help="Use both forward and reverse reads for assembly (if available).")
    parser.add_argument("--skip_preprocessing", action='store_true', help="Skip trimming, host depletion, and Kraken2 on raw reads")
    # Flags for downstream analysis
    parser.add_argument("--bacteria", action='store_true', help="Generate bacterial abundance plots.")
    parser.add_argument("--virus", action='store_true', help="Generate viral abundance plots.")
    parser.add_argument("--archaea", action='store_true', help="Generate archaea abundance plots.")
    parser.add_argument("--eukaryota", action='store_true', help="Generate eukaryota abundance plots.")
    parser.add_argument("--run_ref_base", action="store_true", help="Run reference-based assembly pipeline for each taxon.")
    parser.add_argument("--run_deno_ref", action="store_true", help="Run de novo reference assembly pipeline for each taxon.")
    # New flags for report and file regeneration control
    parser.add_argument("--skip_multiqc", action='store_true', help="Skip MultiQC per user request.")
    parser.add_argument("--skip_reports", action='store_true', help="Skip processing Kraken reports and generating plots; directly run genome assembly.")
    parser.add_argument("--filtered_tsv", help="Path to a filtered merged Kraken output TSV for assembly (optional).")
    parser.add_argument("--skip_existing", action='store_true', help="Skip regenerating files that already exist.")
    parser.add_argument("--process_all_ranks", action='store_true', help="Process all taxonomic ranks.")
    # Filtering arguments
    parser.add_argument("--col_filter", type=str, nargs='+', help="Taxa names to be removed from results.")
    parser.add_argument("--pat_to_keep", type=str, nargs='+', help="Taxa names to specifically keep in results.")
    
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    validate_inputs(args)
    
    # Primary processing: process individual samples.
    process_samples(args)
    
    # Aggregate basic Kraken results using metadata or sample IDs.
    merged_tsv_path = handle_metadata(args)
    
    # ---------------------- Downstream Processing ----------------------
    if not args.skip_multiqc:
        run_multiqc(args.output_dir)
    else:
        logging.info("Skipping MultiQC per user request.")
    
    if args.no_metadata:
        sample_id_df = create_sample_id_df(args.input_dir)
        sample_id_df.to_csv(os.path.join(args.output_dir, "sample_ids.csv"), index=False)
    else:
        sample_id_df = None

    # Define per-domain thresholds (using same values for all domains here for simplicity)
    min_read_counts = {"Bacteria": 1, "Viruses": 1, "Archaea": 1, "Eukaryota": 1}
    max_read_counts = {"Bacteria": args.max_read_count, "Viruses": args.max_read_count, 
                       "Archaea": args.max_read_count, "Eukaryota": args.max_read_count}
    
    run_bowtie = not args.no_bowtie2 and args.bowtie2_index is not None
    
    if not args.skip_reports:
        process_kraken_reports(args.output_dir)
        domains_to_process = []
        if args.bacteria:
            domains_to_process.append("Bacteria")
        if args.virus:
            domains_to_process.append("Viruses")
        if args.archaea:
            domains_to_process.append("Archaea")
        if args.eukaryota:
            domains_to_process.append("Eukaryota")
        
        # Domain rank codes: viruses use a shorter list, others use species, family, domain, kingdom.
        domain_rank_codes = {
            "Bacteria": ['S', 'F', 'D', 'K'],
            "Viruses": ['S', 'F', 'D'],
            "Archaea": ['S', 'F', 'D', 'K'],
            "Eukaryota": ['S', 'F', 'D', 'K']
        }
        
        for domain in domains_to_process:
            logging.info(f"Aggregating results for domain: {domain}")
            for rank in domain_rank_codes.get(domain, ['S']):
                merged_tsv = aggregate_kraken_results(
                    args.output_dir,
                    args.metadata_file,
                    sample_id_df,
                    min_read_counts=min_read_counts,
                    max_read_counts=max_read_counts,
                    rank_code=rank,
                    domain_filter=domain
                )
                # If filtered TSV provided, use it instead.
                if args.filtered_tsv:
                    if os.path.isfile(args.filtered_tsv):
                        merged_tsv = args.filtered_tsv
                        logging.info(f"Using user-provided filtered TSV for assembly: {merged_tsv}")
                    else:
                        logging.error(f"Filtered TSV file '{args.filtered_tsv}' not found.")
                        sys.exit(1)
                if merged_tsv and os.path.isfile(merged_tsv):
                    logging.info(f"Generating abundance plots for {domain} at rank {rank}.")
                    generate_abundance_plots(merged_tsv, args.top_N, args.col_filter, args.pat_to_keep, rank)
                    if args.run_ref_base:
                        df = pd.read_csv(merged_tsv, sep="\t")
                        df = df[~df['Scientific_name'].str.contains('Homo sapiens', case=False, na=False)]
                        df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
                        logging.info(f"Starting reference-based assembly for {domain} at rank {rank}.")
                        ref_based(df, run_bowtie, args.output_dir, skip_existing=args.skip_existing)
                    if args.run_deno_ref:
                        df = pd.read_csv(merged_tsv, sep="\t")
                        df = df[~df['Scientific_name'].str.contains('Homo sapiens', case=False, na=False)]
                        df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
                        logging.info(f"Starting de novo reference assembly for {domain} at rank {rank}.")
                        deno_ref_based(df, args.output_dir, args.output_dir, run_bowtie, skip_existing=args.skip_existing)
    else:
        logging.info("Skipping domain-specific report processing and plot generation. Running genome assembly directly using unfiltered merged TSV.")
        merged_tsv = generate_unfiltered_merged_tsv(args.output_dir, args.metadata_file, sample_id_df)
        if args.filtered_tsv:
            if os.path.isfile(args.filtered_tsv):
                merged_tsv = args.filtered_tsv
                logging.info(f"Using user-provided filtered TSV for assembly: {merged_tsv}")
            else:
                logging.error(f"Filtered TSV file '{args.filtered_tsv}' not found.")
                sys.exit(1)
        if merged_tsv and os.path.isfile(merged_tsv):
            if args.run_ref_base:
                df = pd.read_csv(merged_tsv, sep="\t")
                df = df[~df['Scientific_name'].str.contains('Homo sapiens', case=False, na=False)]
                df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
                logging.info("Starting reference-based assembly for all genomes using provided TSV.")
                ref_based(df, run_bowtie, args.output_dir, skip_existing=args.skip_existing)
            if args.run_deno_ref:
                df = pd.read_csv(merged_tsv, sep="\t")
                df = df[~df['Scientific_name'].str.contains('Homo sapiens', case=False, na=False)]
                df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
                logging.info("Starting de novo reference assembly for all genomes using provided TSV.")
                deno_ref_based(df, args.output_dir, args.output_dir, run_bowtie, skip_existing=args.skip_existing)
    
    # Optionally process all taxonomic ranks
    if args.process_all_ranks:
        if args.no_metadata:
            sample_id_df = create_sample_id_df(args.input_dir)
            sample_id_df.to_csv(os.path.join(args.output_dir, "sample_ids.csv"), index=False)
            process_all_ranks(args.output_dir, sample_id_df=sample_id_df,
                              read_count=args.read_count, max_read_count=args.top_N,
                              top_N=args.top_N, col_filter=args.col_filter, pat_to_keep=args.pat_to_keep)
        else:
            if not args.metadata_file or not os.path.isfile(args.metadata_file):
                logging.error(f"Metadata file '{args.metadata_file}' not found.")
                sys.exit(1)
            process_all_ranks(args.output_dir, metadata_file=args.metadata_file,
                              read_count=args.read_count, max_read_count=args.top_N,
                              top_N=args.top_N, col_filter=args.col_filter, pat_to_keep=args.pat_to_keep)

if __name__ == "__main__":
    main()
