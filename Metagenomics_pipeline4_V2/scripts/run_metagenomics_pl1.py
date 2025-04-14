#!/usr/bin/env python3

import os
import glob
import argparse
import pandas as pd
import sys
import logging
from Metagenomics_pipeline4_V0.kraken_abundance_pipeline import (
    process_sample, 
    aggregate_kraken_results, 
    generate_abundance_plots
)
from Metagenomics_pipeline4_V0.ref_based_assembly import ref_based
from Metagenomics_pipeline4_V0.deno_ref_assembly2 import deno_ref_based

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
        sample_id = sample_id.replace("_R1_001.fastq.gz", "").replace("_R1_001.fastq", "")
        sample_id = sample_id.replace("_R1.fastq.gz", "").replace("_R1.fastq", "")
        sample_id = sample_id.replace("R1.fastq.gz", "").replace("R1.fastq", "")
        sample_id = sample_id.replace("_R1_001", "").replace("_R1", "")  # For cases without ".fastq"
        sample_ids.append(sample_id)

    sample_id_df = pd.DataFrame(sample_ids, columns=["Sample_IDs"])
    return sample_id_df

def validate_inputs(args):
    """Validate input arguments and files"""
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
    """Process all samples according to the specified parameters"""
    run_bowtie = not args.no_bowtie2 and args.bowtie2_index is not None
    
    for forward in glob.glob(os.path.join(args.input_dir, "*_R1*.fastq*")):
        base_name = os.path.basename(forward)
        for pattern in [
            "_R1_001.fastq.gz", "_R1_001.fastq", 
            "_R1.fastq.gz", "_R1.fastq", 
            "R1.fastq.gz", "R1.fastq",
            "_R1_001", "_R1"
        ]:
            base_name = base_name.replace(pattern, "")
        
        # Find matching reverse file if exists
        reverse = None
        if not args.use_assembly or args.paired_assembly:
            reverse_candidates = [
                os.path.join(args.input_dir, f"{base_name}_R2_001.fastq.gz"),
                os.path.join(args.input_dir, f"{base_name}_R2.fastq.gz"),
                os.path.join(args.input_dir, f"{base_name}_R2.fastq"),
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
            skip_preprocessing=args.skip_preprocessing
        )

def handle_metadata(args):
    """Handle metadata processing and result aggregation"""
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
    
    # Flags
    parser.add_argument("--no_bowtie2", action='store_true', help="Skip Bowtie2 host depletion.")
    parser.add_argument("--no_metadata", action='store_true', help="Use sample IDs instead of metadata file.")
    parser.add_argument("--bacteria", action='store_true', help="Generate bacterial abundance plots.")
    parser.add_argument("--virus", action='store_true', help="Generate viral abundance plots.")
    parser.add_argument("--use_precomputed_reports", action='store_true', 
                       help="Use existing Kraken reports instead of running Kraken2.")
    parser.add_argument("--use_assembly", action='store_true', 
                       help="Perform de novo assembly before classification.")
    parser.add_argument("--paired_assembly", action='store_true',
                       help="Use both forward and reverse reads for assembly (if available).")
    parser.add_argument("--run_ref_base", action="store_true", 
                       help="Run reference-based assembly pipeline for each taxon.")
    parser.add_argument("--run_deno_ref", action="store_true", 
                       help="Run denovo reference assembly pipeline for each taxon.")
    parser.add_argument('--skip_preprocessing', action='store_true', help='Skip trimming, host depletion, and Kraken2 on raw reads')
    
    # Filtering arguments
    parser.add_argument("--col_filter", type=str, nargs='+', 
                       help="Taxa names to be removed from results.")
    parser.add_argument("--pat_to_keep", type=str, nargs='+', 
                       help="Taxa names to specifically keep in results.")
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Validate inputs
    validate_inputs(args)
    
    # Process samples
    process_samples(args)
    
    # Handle metadata and aggregate results
    merged_tsv_path = handle_metadata(args)
    
    # Generate abundance plots if requested
    if merged_tsv_path and os.path.isfile(merged_tsv_path):
        if args.virus or args.bacteria:
            logging.info("Generating abundance plots")
            generate_abundance_plots(
                merged_tsv_path, 
                args.top_N,
                args.col_filter, 
                args.pat_to_keep
            )
    
    # Run additional pipelines if requested
    if args.run_ref_base or args.run_deno_ref:
        df = pd.read_csv(merged_tsv_path, sep='\t')
        df = df[df['Scientific_name'].str.contains('virus', case=False, na=False)]
        df = df.apply(lambda col: col.map(lambda x: x.strip() if isinstance(x, str) else x))
        
        if args.run_ref_base:
            logging.info("Starting reference-based pipeline")
            ref_based(df, not args.no_bowtie2, args.output_dir)
            
        if args.run_deno_ref:
            logging.info("Starting denovo reference assembly pipeline")
            deno_ref_based(df, args.output_dir, args.output_dir, not args.no_bowtie2)

if __name__ == "__main__":
    main()
