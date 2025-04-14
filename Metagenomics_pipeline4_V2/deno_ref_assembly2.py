import os
import subprocess
from Bio import SeqIO
import pandas as pd
import logging
import os
from pathlib import Path



# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("pipeline.log"),
        logging.StreamHandler()
    ]
)





def split_fasta(input_file, output_dir):
    """
    Splits a multi-sequence FASTA file into individual FASTA files, 
    each corresponding to a separate sequence, and returns a list of 
    those file paths.

    Args:
    - input_file (str): Path to the input FASTA file.
    - output_dir (str): Directory to store the individual FASTA files.

    Returns:
    - List of paths to the generated FASTA files.
    """
    Path(output_dir).mkdir(exist_ok=True)  # Ensure the output directory exists
    fasta_files = []  # List to store the paths of the generated FASTA files

    with open(input_file, "r") as fasta:
        sequence = []
        accession = None
        for line in fasta:
            if line.startswith(">"):  # New sequence header
                if accession:  # Save the previous sequence
                    file_path = f"{output_dir}/{accession}.fasta"
                    with open(file_path, "w") as out_fasta:
                        out_fasta.write("".join(sequence))
                    fasta_files.append(file_path)  # Add to list of FASTA file paths
                # Extract accession number from header
                accession = line.split()[0][1:]  # Remove '>' and take first part
                sequence = [line]  # Start a new sequence
            else:
                sequence.append(line)
        # Save the last sequence
        if accession:
            file_path = f"{output_dir}/{accession}.fasta"
            with open(file_path, "w") as out_fasta:
                out_fasta.write("".join(sequence))
            fasta_files.append(file_path)  # Add to list of FASTA file paths

    return fasta_files





def get_best_reference(sample_r1, sample_r2, reference_list):
    """
    Align paired-end FASTQ files to a list of reference FASTA files using BWA
    and return the best reference based on alignment scores.

    Parameters:
        sample_r1 (str): Path to the first paired-end FASTQ file.
        sample_r2 (str): Path to the second paired-end FASTQ file.
        reference_list (list): List of paths to reference FASTA files.

    Returns:
        str: The best reference file with the highest alignment score.
    """
    # Dictionary to store alignment scores
    alignment_scores = {}

    for fasta in reference_list:
        index_base = Path(fasta).stem  # Extract the base name for the index
        output_sam = f"{index_base}_aligned.sam"  # SAM file for the output

        # Check if the BWA index exists; if not, create it
        if not Path(f"{fasta}.bwt").exists():
            print(f"Index for {fasta} not found. Creating index...")
            try:
                subprocess.run(["bwa", "index", fasta], check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error in BWA index creation for {fasta}: {e}")
                continue

        # Run BWA for alignment
        bwa_command = [
            "bwa", "mem", fasta, sample_r1, sample_r2
        ]
        try:
            with open(output_sam, "w") as sam_file:
                subprocess.run(bwa_command, check=True, stdout=sam_file)
        except subprocess.CalledProcessError as e:
            print(f"Error in BWA alignment for {fasta}: {e}")
            continue

        # Parse alignment score from the SAM file
        score = 0
        try:
            with open(output_sam, "r") as sam_file:
                for line in sam_file:
                    if not line.startswith("@"):  # Skip header lines
                        fields = line.strip().split("\t")
                        try:
                            # Use the AS:i:<score> tag for alignment score (if available)
                            score += int([tag for tag in fields if tag.startswith("AS:i:")][0].split(":")[2])
                        except IndexError:
                            continue
        except FileNotFoundError:
            print(f"Failed to open SAM file: {output_sam}")
            continue

        alignment_scores[fasta] = score

    # Pick the reference with the highest alignment score
    if alignment_scores:
        best_reference = max(alignment_scores, key=alignment_scores.get)
        print(f"The best reference is {best_reference} with a score of {alignment_scores[best_reference]}")
        return best_reference
    else:
        print("No alignments were successful.")
        return None

def download_and_index_reference(tax, scientific_name, tax_dir):
    """
    Download and index the reference genome for a given taxonomic ID.
    """
    fasta_file = os.path.join(tax_dir, f"{scientific_name}.fasta")
    command = f'esearch -db nucleotide -query "txid{tax}[Organism]" | efilter -source refseq | efetch -format fasta > {fasta_file}'
    try:
        logging.info(f"Downloading reference genome for {scientific_name} (taxid {tax})...")
        subprocess.run(command, shell=True, check=True)
        subprocess.run(f"bwa index {fasta_file}", shell=True, check=True)
        logging.info(f"Reference genome downloaded and indexed: {fasta_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error preparing reference for {tax}: {e}")
        return None
    return fasta_file
def run_de_novo_assembly(sample, sample_r1, sample_r2, output_dir):
    """
    Run de novo assembly using MetaSPAdes.
    """
    contigs_file = os.path.join(f"{output_dir}/{sample}_denovo", "contigs.fasta")
    if os.path.exists(contigs_file):
        logging.info(f"De novo assembly already exists for {sample}. Skipping.")
        return contigs_file

    command = f"metaspades.py -1 {sample_r1} -2 {sample_r2} -o {output_dir}/{sample}_denovo -t 32"
    try:
        logging.info(f"Running de novo assembly for {sample}: {command}")
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in de novo assembly for {sample}: {e}")
        return None

    if not os.path.exists(contigs_file):
        logging.error(f"De novo assembly failed: contigs not generated for {sample}")
        return None

    logging.info(f"De novo assembly completed for {sample}: {contigs_file}")
    return contigs_file

def extract_first_contig_id(fasta_file, output_file):
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                contig_id = line[1:].strip()  # Remove the '>' and any leading/trailing whitespace
                with open(output_file, 'w') as out_file:
                    out_file.write(contig_id + '\n')  # Write contig ID to the output file
                break  # Exit after writing the first contig ID

def calculate_average_read_depth(bam_file):
    """
    Calculate the average read depth from a BAM file using samtools depth.
    """
    try:
        # Run samtools depth command
        result = subprocess.run(
            ["samtools", "depth", bam_file],
            capture_output=True,
            text=True,
            check=True
        )
        # Process the output to calculate average depth
        depths = [int(line.split()[2]) for line in result.stdout.splitlines()]
        if depths:
            average_depth = sum(depths) / len(depths)
        else:
            average_depth = 0
        return average_depth
    except subprocess.CalledProcessError as e:
        logging.error(f"Error calculating read depth for {bam_file}: {e}")
        return None
def deno_ref_based(df, input_dir, output_dir, run_bowtie):
    base_dir = "Fasta_files"

    # Initialize new columns
    df['Ref_len'] = ""
    df['Consensus_len'] = ""
    df['Contig_len'] = ""
    df['Completeness(%)'] = ""
    df['Depth'] = ""
    df['Accession_number'] = ""
    df['sequence'] = ""
    

    dfs = []
    taxids = df['NCBI_ID'].unique()

    for tax in taxids:
        dftax = df[df['NCBI_ID'] == tax].copy()
        scientific_name = dftax['Scientific_name'].iloc[0].replace(' ', '_')
        tax_dir = os.path.join(base_dir, f"{scientific_name}_txid{tax}")
        os.makedirs(tax_dir, exist_ok=True)

        fasta_file = download_and_index_reference(tax, scientific_name, tax_dir)
        if not fasta_file:
            logging.error(f"Skipping {scientific_name} due to missing reference genome.")
            continue

        # Process each sample
        for sample in dftax['SampleID']:
            try:
                sample_r1 = os.path.join(input_dir, f"{sample}_unmapped_1.fastq.gz" if run_bowtie else f"{sample}_trimmed_R1.fastq.gz")
                sample_r2 = os.path.join(input_dir, f"{sample}_unmapped_2.fastq.gz" if run_bowtie else f"{sample}_trimmed_R2.fastq.gz")

                if not (os.path.exists(sample_r1) and os.path.exists(sample_r2)):
                    logging.warning(f"Skipping {sample}: Missing FASTQ files.")
                    continue

                # Run de novo assembly
                contigs_file = run_de_novo_assembly(sample, sample_r1, sample_r2, output_dir)
                if not contigs_file:
                    logging.error(f"Skipping {sample} due to de novo assembly failure.")
                    continue

                # Get best reference from pre-split reference genome
                reference_list = split_fasta(fasta_file, f"{output_dir}/{sample}_{scientific_name}")
                fasta_file = get_best_reference(sample_r1, sample_r2, reference_list)

                # Extract accession number
                cmd = f"grep '^>' {fasta_file} | cut -d ' ' -f1 | sed 's/^>//'"
                acc = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
                acc_id = acc.stdout.strip().split("\n")[0] if acc.stdout.strip() else "Unknown"

                # Run RagTag scaffolding
                rag_dir = os.path.join(output_dir, f"{sample}_{scientific_name}_rag")
                rag_file = os.path.join(rag_dir, "ragtag.scaffold.fasta")
                if not os.path.exists(rag_file):
                    os.makedirs(rag_dir, exist_ok=True)
                    subprocess.run(f"ragtag.py scaffold {fasta_file} {contigs_file} -o {rag_dir}", shell=True, check=True)

                # Extract first contig
                first_contig = os.path.join(rag_dir, "first_contig.fasta")
                output_file = os.path.join(rag_dir, "Output_seqid.txt")
                extract_first_contig_id(rag_file, output_file)
                subprocess.run(f"seqtk subseq {rag_file} {output_file} > {first_contig}", shell=True, check=True)
                subprocess.run(f"bwa index {first_contig}", shell=True, check=True)

                # Mapping & Variant Calling
                sample_dir = os.path.join(base_dir, f"{scientific_name}_assembled")
                os.makedirs(sample_dir, exist_ok=True)
                bam_file = os.path.join(sample_dir, f"{sample}_mapped_reads.bam")
                vcf_file = os.path.join(sample_dir, f"{sample}_variants.vcf")
                consensus_file = os.path.join(sample_dir, f"{sample}_consensus_genome.fa")

                bwa_command = f"bwa mem -a -t 16 {first_contig} {sample_r1} {sample_r2} | samtools view -u - | samtools sort -o {bam_file}"
                subprocess.run(bwa_command, shell=True, check=True)
                subprocess.run(f"samtools index {bam_file}", shell=True, check=True)
                subprocess.run(f"bcftools mpileup -f {first_contig} {bam_file} | bcftools call -c --ploidy 1 -v -o {vcf_file}", shell=True, check=True)
                subprocess.run(f"samtools mpileup -aa -A -d 1000000 -Q 20 -f {first_contig} {bam_file} | ivar consensus -p {consensus_file} -t 0.5 -m 10", shell=True, check=True)

                # Compute sequence lengths and completeness
                def calculate_length(fasta_file):
                    if not os.path.exists(fasta_file):
                        return 0
                    cmd = f"grep -v '^>' {fasta_file} |tr -d '\n' | tr -cd 'ACTG' | wc -c"
                    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
                    return int(result.stdout.strip()) if result.stdout.strip().isdigit() else 0

                ref_len = calculate_length(fasta_file)
                contig_len = calculate_length(first_contig)
                consensus_len = calculate_length(consensus_file)
                completeness = round((consensus_len / ref_len) * 100, 2) if ref_len > 0 else 0
                sequence = str(SeqIO.read(consensus_file, "fasta").seq) if os.path.exists(consensus_file) else ""
                depth=calculate_average_read_depth(bam_file)

                logging.info(f"{sample}: Reference={ref_len}, Contig={contig_len}, Consensus={consensus_len}, Completeness={completeness}%")

                # Update DataFrame
                dftax.loc[dftax['SampleID'] == sample, ['Ref_len', 'Contig_len', 'Consensus_len', 'Completeness(%)','Depth', 'Accession_number', 'sequence']] = \
                    [ref_len, contig_len, consensus_len, completeness,depth, acc_id, sequence]

            except Exception as e:
                logging.error(f"Error processing sample {sample}: {e}")

        dfs.append(dftax)

    # Merge results and save output
    merged_df = pd.concat(dfs, ignore_index=True)
    filtered_df = merged_df[pd.to_numeric(merged_df['Completeness(%)'], errors='coerce') >= 60]
    filtered_df.to_csv("Output-summary_complete.csv", index=False)

    # Remove sequences before saving full summary
    merged_df.drop(columns=['sequence'], inplace=True)
    merged_df.to_csv("Output-summary_all.csv", index=False)
