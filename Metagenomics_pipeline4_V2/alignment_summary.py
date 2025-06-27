import os
import pandas as pd
import subprocess
from Bio import SeqIO

def run_alignment_summary(diamond_tsv: str,
                          merged_fasta: str,
                          fastq_dir: str,
                          output_file: str,
                          tmp_dir: str = "tmp_alignments",
                          run_alignment: bool = True,
                          threads: int = 8,
                          min_contig_len: int = 500) -> str:
    """
    Aligns unmapped reads back to individual contigs using BWA and summarizes mapped read counts per contig.

    Args:
        diamond_tsv (str): Path to diamond_results_contig_with_sampleid.tsv
        merged_fasta (str): Path to merged contigs FASTA
        fastq_dir (str): Directory with unmapped FASTQ files
        output_file (str): Path to write alignment summary TSV
        tmp_dir (str): Temp directory for intermediates
        run_alignment (bool): If False, skip alignments
        threads (int): Threads for BWA
        min_contig_len (int): Min contig length to align

    Returns:
        str: Path to TSV summary
    """
    if not run_alignment:
        print("⚠️ Alignment disabled. Skipping.")
        return ""

    os.makedirs(tmp_dir, exist_ok=True)

    # Load contigs
    contig_dict = {
        rec.id: rec
        for rec in SeqIO.parse(merged_fasta, "fasta")
        if len(rec.seq) >= min_contig_len
    }

    print(f"🔍 Loaded {len(contig_dict)} contigs ≥ {min_contig_len} bp")

    # Load hits
    df = pd.read_csv(diamond_tsv, sep="\t")
    results = []

    for _, row in df.iterrows():
        sample_id = row["Sample_ID"]
        contig_id = row.get("Contig_ID") or row.get("query_id")

        if contig_id not in contig_dict:
            print(f"[!] Contig {contig_id} missing or too short. Skipping.")
            continue

        contig = contig_dict[contig_id]
        safe_id = contig_id.replace("|", "_")
        prefix = f"{sample_id}_{safe_id}"
        ref_fasta = os.path.join(tmp_dir, f"{prefix}.fasta")
        bam_out = os.path.join(tmp_dir, f"{prefix}.bam")

        r1 = os.path.join(fastq_dir, f"{sample_id}_unmapped_1.fastq.gz")
        r2 = os.path.join(fastq_dir, f"{sample_id}_unmapped_2.fastq.gz")
        if not os.path.exists(r1) or not os.path.exists(r2):
            print(f"[!] Missing FASTQ for {sample_id}. Skipping.")
            continue

        # Write reference
        if not os.path.exists(ref_fasta):
            SeqIO.write(contig, ref_fasta, "fasta")

        # Index contig (once)
        if not all(os.path.exists(ref_fasta + ext) for ext in [".bwt", ".amb", ".ann", ".pac", ".sa"]):
            subprocess.run(["bwa", "index", ref_fasta], check=True)

        if os.path.exists(bam_out):
            print(f"📄 Using existing BAM: {bam_out}")
        else:
            try:
                bwa = subprocess.Popen(
                    ["bwa", "mem", "-t", str(threads), ref_fasta, r1, r2],
                    stdout=subprocess.PIPE
                )
                subprocess.run(["samtools", "view", "-bS", "-o", bam_out], stdin=bwa.stdout, check=True)
                bwa.stdout.close()
                bwa.wait()
            except subprocess.CalledProcessError:
                print(f"[!] BWA failed for {sample_id} vs {contig_id}")
                continue

        # Count mapped reads in BAM
        try:
            result = subprocess.run(
                ["samtools", "view", "-c", "-F", "4", bam_out],
                capture_output=True, text=True, check=True
            )
            mapped_reads = int(result.stdout.strip())
        except subprocess.CalledProcessError:
            print(f"[!] BAM read count failed for {prefix}")
            mapped_reads = 0

        results.append({
            "Sample_ID": sample_id,
            "Contig_ID": contig_id,
            "Aligned_Reads": mapped_reads
        })
        print(f"✅ {sample_id} vs {contig_id}: {mapped_reads} reads aligned.")

        # Cleanup
        try:
            os.remove(bam_out)
            os.remove(ref_fasta)
            for ext in [".bwt", ".amb", ".ann", ".pac", ".sa"]:
                index_file = ref_fasta + ext
                if os.path.exists(index_file):
                    os.remove(index_file)
        except Exception as e:
            print(f"[!] Failed to delete intermediary files for {prefix}: {e}")

    # Save summary
    summary_df = pd.DataFrame(results)
    summary_df.to_csv(output_file, sep="\t", index=False)
    print(f"\n📝 Alignment summary saved to: {output_file}")

    return output_file
