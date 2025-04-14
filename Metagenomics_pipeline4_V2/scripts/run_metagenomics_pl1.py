
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
    parser.add_argument("--kraken_db", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--input_dir", required=True)
    parser.add_argument("--bowtie2_index")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--metadata_file")
    parser.add_argument("--read_count", type=int, default=1)
    parser.add_argument("--top_N", type=int, default=10000)
    parser.add_argument("--max_read_count", type=int, default=5000000000)
    parser.add_argument("--no_bowtie2", action='store_true')
    parser.add_argument("--no_metadata", action='store_true')
    parser.add_argument("--use_precomputed_reports", action='store_true')
    parser.add_argument("--use_assembly", action='store_true')
    parser.add_argument("--paired_assembly", action='store_true')
    parser.add_argument("--skip_preprocessing", action='store_true')
    parser.add_argument("--bacteria", action='store_true')
    parser.add_argument("--virus", action='store_true')
    parser.add_argument("--archaea", action='store_true')
    parser.add_argument("--eukaryota", action='store_true')
    parser.add_argument("--run_ref_base", action="store_true")
    parser.add_argument("--run_deno_ref", action="store_true")
    parser.add_argument("--skip_multiqc", action='store_true')
    parser.add_argument("--skip_reports", action='store_true')
    parser.add_argument("--filtered_tsv")
    parser.add_argument("--skip_existing", action='store_true')
    parser.add_argument("--process_all_ranks", action='store_true')
    parser.add_argument("--col_filter", type=str, nargs='+')
    parser.add_argument("--pat_to_keep", type=str, nargs='+')

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    validate_inputs(args)

    process_samples(args)
    merged_tsv_path = handle_metadata(args)

    if not args.skip_multiqc:
        run_multiqc(args.output_dir)

    if args.no_metadata:
        sample_id_df = create_sample_id_df(args.input_dir)
    else:
        sample_id_df = None

    if not args.skip_reports:
        process_kraken_reports(args.output_dir)
        domains = ["Bacteria", "Viruses", "Archaea", "Eukaryota"]
        domain_flags = [args.bacteria, args.virus, args.archaea, args.eukaryota]
        domain_rank_codes = {"Bacteria": ['S', 'F', 'D', 'K'], "Viruses": ['S', 'F', 'D'], "Archaea": ['S', 'F', 'D', 'K'], "Eukaryota": ['S', 'F', 'D', 'K']}

        for domain, flag in zip(domains, domain_flags):
            if flag:
                for rank in domain_rank_codes[domain]:
                    merged_tsv = aggregate_kraken_results(
                        args.output_dir, args.metadata_file, sample_id_df,
                        {domain: args.read_count}, {domain: args.max_read_count}, rank, domain
                    )
                    if args.filtered_tsv and os.path.isfile(args.filtered_tsv):
                        merged_tsv = args.filtered_tsv
                    if merged_tsv and os.path.isfile(merged_tsv):
                        generate_abundance_plots(merged_tsv, args.top_N, args.col_filter, args.pat_to_keep, rank)
                        if args.run_ref_base:
                            df = pd.read_csv(merged_tsv, sep="\t")
                            df = df[~df['Scientific_name'].str.contains('Homo sapiens', case=False, na=False)]
                            df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
                            ref_based(df, args.output_dir, args.input_dir, args.bowtie2_index, args.threads, rank)
                        if args.run_deno_ref:
                            deno_ref_based(merged_tsv, args.output_dir, args.input_dir, args.threads, rank)
    elif args.run_ref_base or args.run_deno_ref:
        merged_tsv = args.filtered_tsv if args.filtered_tsv else merged_tsv_path
        if os.path.isfile(merged_tsv):
            if args.run_ref_base:
                ref_based(pd.read_csv(merged_tsv, sep="\t"), args.output_dir, args.input_dir, args.bowtie2_index, args.threads, "S")
            if args.run_deno_ref:
                deno_ref_based(merged_tsv, args.output_dir, args.input_dir, args.threads, "S")

if __name__ == "__main__":
    main()
