#!/usr/bin/env python3
"""
run_metagenomics_pl1.py

Metagenomics pipeline including Kraken2 classification, optional assembly‐based contig
extraction & Diamond annotation, and downstream analyses.
"""

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
    process_all_ranks,
    process_output_reports
)
from Metagenomics_pipeline4_V2.ref_based_assembly import ref_based
from Metagenomics_pipeline4_V2.deno_ref_assembly2 import deno_ref_based

from Metagenomics_pipeline4_V2 import extract_contigs_diamond  # our contig‐extraction & Diamond module

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
        sid = os.path.basename(f)
        for pat in ["_R1_001.fastq.gz", "_R1_001.fastq", "_R1.fastq.gz",
                    "_R1.fastq", "R1.fastq.gz", "R1.fastq", "_R1_001", "_R1"]:
            sid = sid.replace(pat, "")
        sample_ids.append(sid)
    return pd.DataFrame(sample_ids, columns=["Sample_IDs"])

def validate_inputs(args):
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

def process_samples(args):
    run_bowtie = not args.no_bowtie2 and args.bowtie2_index is not None
    for forward in glob.glob(os.path.join(args.input_dir, "*_R1*.fastq*")):
        base = os.path.basename(forward)
        for pat in ["_R1_001.fastq.gz", "_R1_001.fastq", "_R1.fastq.gz",
                    "_R1.fastq", "R1.fastq.gz", "R1.fastq", "_R1_001", "_R1"]:
            base = base.replace(pat, "")
        reverse = None
        if not args.use_assembly or args.paired_assembly:
            candidates = [
                os.path.join(args.input_dir, f"{base}_R2_001.fastq.gz"),
                os.path.join(args.input_dir, f"{base}_R2.fastq.gz"),
                os.path.join(args.input_dir, f"{base}_R2.fastq")
            ]
            reverse = next((r for r in candidates if os.path.isfile(r)), None)
            if not reverse and not args.use_assembly:
                logging.warning(f"No R2 found for {base}, skipping.")
                continue

        logging.info(f"Processing sample {base} (assembly={args.use_assembly})")
        process_sample(
            forward=forward,
            reverse=reverse,
            base_name=base,
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

def main():
    parser = argparse.ArgumentParser(
        description="Metagenomics pipeline for taxonomic classification and analysis"
    )
    parser.add_argument("--kraken_db", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--input_dir", required=True)
    parser.add_argument("--bowtie2_index")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument(
        "--diamond_db",
        default="/home/soumareh/mnt/nrdb/nr",
        help="Path to the Diamond-formatted NR database"
    )
    parser.add_argument("--metadata_file")
    parser.add_argument("--read_count", type=int, default=1)
    parser.add_argument("--top_N", type=int, default=10000)
    parser.add_argument("--max_read_count", type=int, default=5000000000)
    parser.add_argument("--no_bowtie2", action="store_true")
    parser.add_argument("--no_metadata", action="store_true")
    parser.add_argument("--use_precomputed_reports", action="store_true")
    parser.add_argument("--use_assembly", action="store_true")
    parser.add_argument("--paired_assembly", action="store_true")
    parser.add_argument("--skip_preprocessing", action="store_true")
    parser.add_argument("--bacteria", action="store_true")
    parser.add_argument("--virus", action="store_true")
    parser.add_argument("--archaea", action="store_true")
    parser.add_argument("--eukaryota", action="store_true")
    parser.add_argument("--run_ref_base", action="store_true")
    parser.add_argument("--run_deno_ref", action="store_true")
    parser.add_argument("--skip_multiqc", action="store_true")
    parser.add_argument("--skip_reports", action="store_true")
    parser.add_argument("--filtered_tsv")
    parser.add_argument("--skip_existing", action="store_true")
    parser.add_argument("--process_all_ranks", action="store_true")
    parser.add_argument("--col_filter", type=str, nargs="+")
    parser.add_argument("--pat_to_keep", type=str, nargs="+")

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    validate_inputs(args)

    process_samples(args)
    merged_tsv = handle_metadata(args)

    if not args.skip_reports:
        logging.info("Processing Kraken reports…")
        process_kraken_reports(args.output_dir)
        logging.info("Processing output reports…")
        process_output_reports(args.output_dir)

    # Assembly‐based contig extraction & Diamond annotation
    if args.use_assembly:
        logging.info("Running contig extraction and Diamond annotation…")
        extract_contigs_diamond.extract_contigs(
            base_contigs_dir=args.output_dir,
            summary_filename=os.path.join(args.output_dir, "contigs_summary.tsv")
        )
        extract_contigs_diamond.merge_and_rename_contigs(
            base_contigs_dir=args.output_dir,
            merged_filename=os.path.join(args.output_dir, "merged_contigs_renamed.fasta")
        )
        extract_contigs_diamond.run_diamond(
            diamond_db=args.diamond_db,
            query_file=os.path.join(args.output_dir, "merged_contigs_renamed.fasta"),
            output_file=os.path.join(args.output_dir, "results.m8"),
            threads=args.threads
        )
        extract_contigs_diamond.process_diamond_results(
            results_filename=os.path.join(args.output_dir, "results.m8"),
            extracted_csv=os.path.join(args.output_dir, "extracted_virus.csv"),
            extracted_csv1=os.path.join(args.output_dir, "extracted_virus1.csv")
        )
        logging.info("Contig extraction and Diamond steps complete.")

    if not args.skip_multiqc:
        run_multiqc(args.output_dir)

    if args.no_metadata:
        sample_id_df = create_sample_id_df(args.input_dir)
    else:
        sample_id_df = None

    if not args.skip_reports:
        domains = ["Bacteria", "Viruses", "Archaea", "Eukaryota"]
        domain_flags = [args.bacteria, args.virus, args.archaea, args.eukaryota]
        domain_rank_codes = {"Bacteria": ['S','S1','S1','S2' ,'F','F1','F2','F3', 'D','D1','D2','D3'], "Viruses": ['S','S1','S1','S2' ,'F','F1','F2','F3', 'D','D1','D2','D3'], "Archaea": ['S','S1','S1','S2' ,'F','F1','F2','F3', 'D','D1','D2','D3'], "Eukaryota": ['S','S1','S1','S2' ,'F','F1','F2','F3', 'D','D1','D2','D3']}

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
