import os
import pandas as pd
import subprocess
from Bio import SeqIO

def run_alignment_summary(diamond_tsv: str,
                          merged_fasta: str,
                          fastq_dir: str,
                          output_file: str,
                          tmp_dir: str = "tmp_alignments",
                          run_alignment: bool = True) -> str:
    """
    Aligns unmapped reads back to contigs using BWA and summarizes mapped read counts.

    Args:
        diamond_tsv (str): Path to diamond_results_contig_with_sampleid.tsv
        merged_fasta (str): Path to the merged viral contigs FASTA file
        fastq_dir (str): Directory containing unmapped FASTQ files
        output_file (str): Path to write alignment summary TSV
        tmp_dir (str): Directory to store intermediate files
        run_alignment (bool): Whether to run the alignment step or skip it

    Returns:
        str: Path to the final summary TSV file
    """
    if not run_alignment:
        print("⚠️ Alignment flag is disabled. Skipping alignment step.")
        return ""

    os.makedirs(tmp_dir, exist_ok=True)

    # Load Diamond hits
    df = pd.read_csv(diamond_tsv, sep="\t")

    # Load all contigs from merged FASTA into a dictionary
    print("🔍 Loading contigs...")
    contig_dict = SeqIO.to_dict(SeqIO.parse(merged_fasta, "fasta"))

    results = []

    for _, row in df.iterrows():
        sample_id = row["Sample_ID"]
        contig_id = row["Contig_ID"] if "Contig_ID" in row else row.get("query_id")

        if contig_id not in contig_dict:
            print(f"[!] Contig {contig_id} not found in FASTA. Skipping.")
            continue

        safe_contig_id = contig_id.replace("|", "_")
        prefix = f"{sample_id}_{safe_contig_id}"
        ref_fasta = os.path.join(tmp_dir, f"{prefix}.fasta")
        sam_out = os.path.join(tmp_dir, f"{prefix}.sam")

        # Write contig to reference FASTA
        SeqIO.write(contig_dict[contig_id], ref_fasta, "fasta")

        # Index with BWA
        subprocess.run(["bwa", "index", ref_fasta], check=True)

        r1 = os.path.join(fastq_dir, f"{sample_id}_unmapped_1.fastq.gz")
        r2 = os.path.join(fastq_dir, f"{sample_id}_unmapped_2.fastq.gz")
        if not os.path.exists(r1) or not os.path.exists(r2):
            print(f"[!] FASTQ files missing for {sample_id}. Skipping.")
            continue

        try:
            with open(sam_out, "w") as samfile:
                subprocess.run(["bwa", "mem", ref_fasta, r1, r2], stdout=samfile, check=True)
        except subprocess.CalledProcessError:
            print(f"[!] BWA alignment failed for {sample_id} vs {contig_id}.")
            continue

        mapped_reads = 0
        with open(sam_out, "r") as f:
            for line in f:
                if line.startswith("@"): continue
                fields = line.split("\t")
                flag = int(fields[1])
                if not flag & 0x4:
                    mapped_reads += 1

        results.append({
            "Sample_ID": sample_id,
            "Contig_ID": contig_id,
            "Aligned_Reads": mapped_reads
        })
        print(f"✅ {sample_id} vs {contig_id} → {mapped_reads} reads aligned.")

    summary_df = pd.DataFrame(results)
    summary_df.to_csv(output_file, sep="\t", index=False)
    print(f"\n📝 Alignment summary saved to: {output_file}")

    return output_file
