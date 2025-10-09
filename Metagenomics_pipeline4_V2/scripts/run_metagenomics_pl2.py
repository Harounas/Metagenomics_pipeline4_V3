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
from Metagenomics_pipeline4_V2.virus_scaffolding import scaffold_virus_contigs

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
        for pat in [
            "_R1_001.fastq.gz", "_R1_001.fastq", "_R1.fastq.gz",
            "_R1.fastq", "R1.fastq.gz", "R1.fastq", "_R1_001", "_R1"
        ]:
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
   
    if args.bowtie2_index:
      bt2 = args.bowtie2_index + ".1.bt2"
      bt2l = args.bowtie2_index + ".1.bt2l"
      if not (os.path.exists(bt2) or os.path.exists(bt2l)):
        logging.error(f"Bowtie2 index '{args.bowtie2_index}' not found (.bt2 or .bt2l).")
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

def process_samples(args):
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

# Truncated main() section and full Diamond + geNomad + scaffolding + alignment logic due to space.
# Will continue in a follow-up script below.

#parser = argparse.ArgumentParser(description="Metagenomics pipeline for taxonomic classification and analysis")
#args = parser.parse_args()

def main():
    parser = argparse.ArgumentParser(description="Metagenomics pipeline for taxonomic classification and analysis")
    # [Arguments truncated for brevity — same as part 1, must be consistent]
     #parser = argparse.ArgumentParser(description="Metagenomics pipeline for taxonomic classification and analysis")

    # Input/Output
    parser.add_argument("--input_dir", required=True, help="Directory with input FASTQ files")
    parser.add_argument("--output_dir", required=True, help="Directory to store outputs")
    parser.add_argument("--metadata_file", help="Metadata file with sample information")
    parser.add_argument("--filtered_tsv", help="Optional filtered TSV file for report and downstream processing")

    # Databases
    parser.add_argument("--kraken_db", required=True, help="Path to Kraken2 database")
    parser.add_argument("--bowtie2_index", help="Path prefix to Bowtie2 host index")
    parser.add_argument("--genomad_db", help="Path to geNomad database")
    parser.add_argument("--diamond_db", help="Path to DIAMOND database (.dmnd)")
    parser.add_argument("--nr_path", help="Path to NR protein FASTA for annotation")

    # Pipeline Controls
    parser.add_argument("--skip_preprocessing", action="store_true", help="Skip fastp trimming")
    parser.add_argument("--no_bowtie2", action="store_true", help="Skip host depletion with Bowtie2")
    parser.add_argument("--use_precomputed_reports", action="store_true", help="Use existing Kraken2 reports")
    parser.add_argument("--use_assembly", action="store_true", help="Enable SPAdes assembly + contig processing")
    parser.add_argument("--skip_existing", action="store_true", help="Skip samples with existing output")
    parser.add_argument("--skip_reports", action="store_true", help="Skip Kraken report generation and plotting")
    parser.add_argument("--skip_diamond", action="store_true", help="Skip Diamond annotation")
    parser.add_argument("--skip_genomad", action="store_true", help="Skip geNomad detection")
    parser.add_argument("--run_genomad", action="store_true", help="Run geNomad on viral contigs")
    parser.add_argument("--run_scaffolding", action="store_true", help="Run scaffolding of viral contigs")
    parser.add_argument("--run_alignment", action="store_true", help="Run alignment summary on viral contigs")
    parser.add_argument("--run_ref_base", action="store_true", help="Run reference-based genome assembly")
    parser.add_argument("--run_deno_ref", action="store_true", help="Run de novo-assisted reference assembly")
    parser.add_argument("--skip_multiqc", action="store_true", help="Skip MultiQC summary report")

    # Domain-specific Kraken plots
    parser.add_argument("--bacteria", action="store_true", help="Generate plots for Bacteria")
    parser.add_argument("--virus", action="store_true", help="Generate plots for Viruses")
    parser.add_argument("--archaea", action="store_true", help="Generate plots for Archaea")
    parser.add_argument("--eukaryota", action="store_true", help="Generate plots for Eukaryota")

    # Plotting/Filtering
    parser.add_argument("--top_N", type=int, default=15, help="Number of top taxa to plot")
    parser.add_argument("--col_filter", help="Column filter for abundance plots")
    parser.add_argument("--pat_to_keep", help="Regex pattern for taxa to keep in plots")

    # Sample and Metadata Handling
    parser.add_argument("--no_metadata", action="store_true", help="No metadata; use filenames as sample IDs")
    parser.add_argument("--read_count", type=int, default=10000, help="Min read count for inclusion")
    parser.add_argument("--max_read_count", type=int, default=10000000, help="Max read count for inclusion")

    # Threading/Concurrency
    parser.add_argument("--threads", type=int, default=8, help="Threads for Kraken, Bowtie2, Diamond, etc.")
    parser.add_argument("--parallel", type=int, default=4, help="Number of samples to process in parallel")
    parser.add_argument("--max_assemblies", type=int, default=2, help="Max concurrent SPAdes assemblies")
    parser.add_argument("--max_workers", type=int, default=4, help="Max worker processes for alignment summary")
    parser.add_argument("--bwa_threads", type=int, default=4, help="Threads per alignment job")

    # Diamond toggle
    parser.add_argument("--diamond", action="store_true", help="Enable Diamond annotation")

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

    if args.diamond and not args.skip_diamond:
        if not args.nr_path:
            logging.error("Missing --nr_path required for Diamond annotation.")
            sys.exit(1)

        if args.use_assembly:
            logging.info("Running contig extraction and Diamond annotation…")
            long_contigs_tsv = os.path.join(args.output_dir, "long_contigs_summary.tsv")
            long_contigs_fasta = os.path.join(args.output_dir, "long_contigs.fasta")
            extract_contigs_diamond.extract_long_contigs_kraken(
                base_contigs_dir=args.output_dir,
                output_tsv=long_contigs_tsv
            )
            long_contigs_records = []
            with open(long_contigs_tsv) as tsv_file:
                reader = csv.DictReader(tsv_file, delimiter="\t")
                for row in reader:
                    sample_id, gene_id = row["Sample_ID"], row["gene"]
                    contig_path = os.path.join(args.output_dir, sample_id, "contigs.fasta")
                    if os.path.exists(contig_path):
                        for rec in SeqIO.parse(contig_path, "fasta"):
                            if rec.id == gene_id:
                                rec.id = f"{sample_id}|{gene_id}"
                                rec.description = ""
                                long_contigs_records.append(rec)
                                break
            SeqIO.write(long_contigs_records, long_contigs_fasta, "fasta")

            if args.run_genomad and not args.skip_genomad:
                genomad_input_fasta = os.path.join(args.output_dir, "merged_contigs_genomad.fasta")
                genomad_out_dir = os.path.join(args.output_dir, "genomad_output")
                extract_contigs_diamond.extract_and_merge_contigs_genomad(
                    base_contigs_dir=args.output_dir,
                    output_fasta=genomad_input_fasta
                )
                genomad_output_viral_fasta = extract_contigs_diamond.run_genomad(
                    input_fasta=genomad_input_fasta,
                    output_dir=genomad_out_dir,
                    genomad_db=args.genomad_db,
                    threads=args.threads
                )
            else:
                logging.warning("Skipping geNomad run; using only Kraken long contigs")
                genomad_output_viral_fasta = long_contigs_fasta

            combined_fasta_for_clustering = os.path.join(args.output_dir, "combined_contigs_for_clustering.fasta")
            extract_contigs_diamond.filter_and_merge(
                fasta_paths=[genomad_output_viral_fasta, long_contigs_fasta],
                min_length=200,
                output_path=combined_fasta_for_clustering
            )

            clustered_out_dir = os.path.join(args.output_dir, "clustered_output")
            clustered_fasta = extract_contigs_diamond.cluster_contigs(
                virus_fasta=combined_fasta_for_clustering,
                output_dir=clustered_out_dir,
                threads=args.threads
            )
            final_long_clustered_fasta = os.path.join(args.output_dir, "clustered_long_contigs.fasta")
            extract_contigs_diamond.extract_long_contigs(
                input_fasta=clustered_fasta,
                output_fasta=final_long_clustered_fasta
            )

            diamond_result_file = os.path.join(args.output_dir, "results_clustered.m8")
            extract_contigs_diamond.run_diamond(
                diamond_db=args.diamond_db,
                query_file=final_long_clustered_fasta,
                output_file=diamond_result_file,
                threads=args.threads
            )

            processed_output = process_virus_contigs(
                fasta_file=args.nr_path,
                diamond_results_file=diamond_result_file,
                output_dir=args.output_dir
            )

            filtered_clusters_file = process_clustered_contigs(
                clstr_file=os.path.join(args.output_dir, "clustered_output", "clustered_contigs.fasta.clstr"),
                diamond_tsv=os.path.join(args.output_dir, "diamond_results_contig_with_sampleid.tsv"),
                output_dir=args.output_dir
            )

            if args.run_scaffolding:
                df = pd.read_csv(os.path.join(args.output_dir, "filtered_clusters_assigned_rep_virus.tsv"), sep="\t")
                unique_pairs = df[["Sample_ID", "virus"]].drop_duplicates()
                for _, row in unique_pairs.iterrows():
                    sample_id = row["Sample_ID"]
                    virus = row["virus"]
                    try:
                        fasta_out, length_out = scaffold_virus_contigs(
                            tsv_path=filtered_clusters_file,
                            sample_id=sample_id,
                            virus_name=virus,
                            contigs_root=args.output_dir,
                            output_root=os.path.join(args.output_dir, "scaffolded_out"),
                            threads=args.threads
                        )
                        logging.info(f"✅ Scaffolded {sample_id} – {virus}: {fasta_out}")
                    except Exception as e:
                        logging.error(f"❌ Scaffold failed for {sample_id} – {virus}: {e}")

            if args.run_alignment:
                logging.info("🧬 Running alignment summary for viral contigs...")
                run_alignment_summary(
                    diamond_tsv=filtered_clusters_file,
                    merged_fasta=combined_fasta_for_clustering,
                    fastq_dir=args.output_dir,
                    output_file=os.path.join(args.output_dir, "alignment_summary.tsv"),
                    tmp_dir=os.path.join(args.output_dir, "tmp_alignments"),
                    run_alignment=args.run_alignment,
                    max_workers=args.max_workers,
                    bwa_threads_per_job=args.bwa_threads,
                    min_contig_len=200
                )
        else:
            logging.warning("Diamond requested but --use_assembly not provided. Skipping contig analysis.")
    else:
        logging.info("⚠️ Skipping Diamond step as requested.")

    
    if not args.skip_reports:
        domains = ["Bacteria", "Viruses", "Archaea", "Eukaryota"]
        domain_flags = [args.bacteria, args.virus, args.archaea, args.eukaryota]
        domain_rank_codes = {
            "Bacteria":  ['S','S1','S2','F','F1','F2','F3','D','D1','D2','D3'],
            "Viruses":   ['S','S1','S2','F','F1','F2','F3','D','D1','D2','D3'],
            "Archaea":   ['S','S1','S2','F','F1','F2','F3','D','D1','D2','D3'],
            "Eukaryota": ['S','S1','S2','F','F1','F2','F3','D','D1','D2','D3']
        }

        sample_id_df = create_sample_id_df(args.input_dir) if args.no_metadata else None

        for domain, flag in zip(domains, domain_flags):
            if flag:
                for rank in domain_rank_codes[domain]:
                    merged_tsv = aggregate_kraken_results(
                        args.output_dir,
                        args.metadata_file,
                        sample_id_df,
                        {domain: args.read_count},
                        {domain: args.max_read_count},
                        rank,
                        domain
                    )
                    if args.filtered_tsv and os.path.isfile(args.filtered_tsv):
                        merged_tsv = args.filtered_tsv
                    if merged_tsv and os.path.isfile(merged_tsv):
                        generate_abundance_plots(merged_tsv, args.top_N, args.col_filter, args.pat_to_keep, rank)
                        if args.run_ref_base:
                            df = pd.read_csv(merged_tsv, sep="\t")
                            df = df[~df['Scientific_name'].str.contains('Homo sapiens', na=False)]
                            df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
                            ref_based(df, args.output_dir, args.input_dir, args.bowtie2_index, args.threads, rank)
                        if args.run_deno_ref:
                            deno_ref_based(merged_tsv, args.output_dir, args.input_dir, args.threads, rank)
    elif args.run_ref_base or args.run_deno_ref:
        merged_tsv_file = args.filtered_tsv if args.filtered_tsv else merged_tsv
        if os.path.isfile(merged_tsv_file):
            if args.run_ref_base:
                ref_based(pd.read_csv(merged_tsv_file, sep="\t"), args.output_dir, args.input_dir, args.bowtie2_index, args.threads, "S")
            if args.run_deno_ref:
                deno_ref_based(merged_tsv_file, args.output_dir, args.input_dir, args.threads, "S")


#if not args.skip_multiqc:
        #run_multiqc(args.output_dir)

if __name__ == "__main__":
    main()
