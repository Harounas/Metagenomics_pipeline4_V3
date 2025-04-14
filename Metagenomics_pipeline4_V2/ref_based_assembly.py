import os
import subprocess
import pandas as pd
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

REFERENCE_DIR = "Reference_FASTA"

def ensure_directory_exists(directory):
    """Creates a directory if it does not exist."""
    Path(directory).mkdir(parents=True, exist_ok=True)

def split_fasta(input_file, output_dir):
    """Splits a multi-sequence FASTA file into individual FASTA files."""
    ensure_directory_exists(output_dir)
    fasta_files = []
    with open(input_file, "r") as fasta:
        sequence = []
        accession = None
        for line in fasta:
            if line.startswith(">"):
                if accession:
                    file_path = f"{output_dir}/{accession}.fasta"
                    with open(file_path, "w") as out_fasta:
                        out_fasta.write("".join(sequence))
                    fasta_files.append(file_path)
                accession = line.split()[0][1:]
                sequence = [line]
            else:
                sequence.append(line)
        if accession:
            file_path = f"{output_dir}/{accession}.fasta"
            with open(file_path, "w") as out_fasta:
                out_fasta.write("".join(sequence))
            fasta_files.append(file_path)
    return fasta_files

def calculate_average_read_depth(bam_file):
    """Calculate average read depth from a BAM file using samtools depth."""
    try:
        result = subprocess.run(["samtools", "depth", bam_file], capture_output=True, text=True, check=True)
        depths = [int(line.split()[2]) for line in result.stdout.splitlines()]
        return sum(depths) / len(depths) if depths else 0
    except subprocess.CalledProcessError as e:
        logging.error(f"Error calculating read depth for {bam_file}: {e}")
        return None

def get_best_reference(sample_r1, sample_r2, reference_list):
    """Aligns paired-end FASTQ files to a list of reference FASTA files and returns the best reference."""
    alignment_scores = {}
    for fasta in reference_list:
        index_base = Path(fasta).stem
        output_sam = f"{index_base}_aligned.sam"

        if not Path(f"{fasta}.bwt").exists():
            subprocess.run(["bwa", "index", fasta], check=True)

        try:
            with open(output_sam, "w") as sam_file:
                subprocess.run(["bwa", "mem", fasta, sample_r1, sample_r2], check=True, stdout=sam_file)
        except subprocess.CalledProcessError as e:
            logging.error(f"Error in BWA alignment for {fasta}: {e}")
            continue

        score = 0
        try:
            with open(output_sam, "r") as sam_file:
                for line in sam_file:
                    if not line.startswith("@"):
                        fields = line.strip().split("\t")
                        try:
                            score += int([tag for tag in fields if tag.startswith("AS:i:")][0].split(":")[2])
                        except IndexError:
                            continue
        except FileNotFoundError:
            logging.error(f"Failed to open SAM file: {output_sam}")
            continue

        alignment_scores[fasta] = score

    return max(alignment_scores, key=alignment_scores.get) if alignment_scores else None

def extract_sequence(fasta_file):
    """Extracts sequence from a FASTA file."""
    try:
        with open(fasta_file, "r") as f:
            return "".join(line.strip() for line in f if not line.startswith(">"))
    except Exception as e:
        logging.error(f"Error extracting sequence from {fasta_file}: {e}")
        return ""

def fetch_reference_fasta(taxid, output_path):
    """Fetches reference FASTA from NCBI if not already downloaded."""
    if os.path.exists(output_path):
        logging.info(f"Reference FASTA for taxid {taxid} already exists.")
        return output_path

    ensure_directory_exists(REFERENCE_DIR)
    command = f'esearch -db nucleotide -query "txid{taxid}[Organism]" | efilter -source refseq | efetch -format fasta > {output_path}'
    subprocess.run(command, shell=True, check=True)
    
    if not os.path.exists(output_path):
        logging.error(f"FASTA file {output_path} was not created.")
        return None

    subprocess.run(f"bwa index {output_path}", shell=True, check=True)
    return output_path

def ref_based(df, run_bowtie, input_dir):
    """Performs reference-based analysis for each taxon in the dataset."""
    ensure_directory_exists(REFERENCE_DIR)
    df['Ref_len'] = ""
    df['Consensus_len'] = ""
    df['Completeness(%)'] = ""
    df['Depth'] = ""
    df['Accession_number'] = ""
    df['sequence'] = ""
    
    results = []
    
    for tax in df['NCBI_ID'].unique():
        dftax = df[df['NCBI_ID'] == tax].copy()
        scientific_name = dftax['Scientific_name'].iloc[0].replace(' ', '_').replace('/', '_').replace(')', '').replace('(', '')
        fasta_path = os.path.join(REFERENCE_DIR, f"{scientific_name}_txid{tax}.fasta")
        
        fasta_file = fetch_reference_fasta(tax, fasta_path)
        if not fasta_file:
            continue
        
        for sample in dftax['SampleID']:
            sample_r1 = os.path.join(input_dir, f"{sample}_unmapped_1.fastq.gz" if run_bowtie else f"{sample}_trimmed_R1.fastq.gz")
            sample_r2 = os.path.join(input_dir, f"{sample}_unmapped_2.fastq.gz" if run_bowtie else f"{sample}_trimmed_R2.fastq.gz")
            
            sample_dir = os.path.join("Fasta_files", f"{scientific_name}_assembled")
            ensure_directory_exists(sample_dir)
            
            reference_list = split_fasta(fasta_file, sample_dir)
            best_ref = get_best_reference(sample_r1, sample_r2, reference_list)
            if not best_ref:
                continue
            
            acc_cmd = f"grep '^>' {best_ref} | cut -d ' ' -f1 | sed 's/^>//'"
            acc = subprocess.run(acc_cmd, shell=True, capture_output=True, text=True, check=True).stdout.strip().split("\n")[0]
            
            bam_file = os.path.join(sample_dir, f"{sample}_{scientific_name}_mapped_reads.bam")
            vcf_file = os.path.join(sample_dir, f"{sample}_{scientific_name}_variants.vcf")
            consensus_file = os.path.join(sample_dir, f"{sample}_{scientific_name}_consensus_genome.fa")
            
            subprocess.run(f"bwa mem -a -t 16 {best_ref} {sample_r1} {sample_r2} | samtools view -u -@ 3 - | samtools sort -@ 16 -o {bam_file}", shell=True, check=True)
            subprocess.run(f"samtools index {bam_file}", shell=True, check=True)
            subprocess.run(f"bcftools mpileup -f {best_ref} {bam_file} | bcftools call -c --ploidy 1 -v -o {vcf_file}", shell=True, check=True)
            subprocess.run(f"samtools mpileup -aa -A -d 0 -Q 0 -f {best_ref} {bam_file} | ivar consensus -p {consensus_file} -t 0.5 -m 10 -n N", shell=True, check=True)

            ref_len = len(extract_sequence(best_ref))
            consensus_len = len(extract_sequence(consensus_file))
            completeness = round((consensus_len / ref_len) * 100, 2) if ref_len > 0 else 0
            depth = calculate_average_read_depth(bam_file)
            
            dftax.loc[dftax['SampleID'] == sample, ['Ref_len', 'Consensus_len', 'Completeness(%)', 'Depth', 'Accession_number', 'sequence']] = [ref_len, consensus_len, completeness, depth, acc, extract_sequence(consensus_file)]
        
        results.append(dftax)
    
    pd.concat(results).to_csv("Output-summary.csv", index=False)

    merged_df = pd.concat(results, ignore_index=True)
    
    
    filtered_df = merged_df[pd.to_numeric(merged_df['Completeness(%)'], errors='coerce') >= 60]
    filtered_df.to_csv("Output-summary_complete.csv", index=False)
    merged_df.drop(columns=['sequence'], inplace=True)
    merged_df.to_csv("Output-summary1.csv", index=False)
  
