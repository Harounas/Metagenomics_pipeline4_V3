import os
import pandas as pd
import subprocess
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, Any, Optional, Tuple

def run_alignment_summary(diamond_tsv: str,
                          merged_fasta: str,
                          fastq_dir: str,
                          output_file: str,
                          tmp_dir: str = "tmp_alignments",
                          run_alignment: bool = True,
                          max_workers: Optional[int] = None,
                          bwa_threads_per_job: int = 4,
                          min_contig_len: int = 500) -> str:
    """
    Aligns unmapped reads back to individual contigs using BWA and summarizes mapped read counts,
    breadth coverage, and mean depth per contig. Runs jobs in parallel.

    Args:
        diamond_tsv (str): Path to diamond_results_contig_with_sampleid.tsv
        merged_fasta (str): Path to merged contigs FASTA
        fastq_dir (str): Directory with unmapped FASTQ files (expects {Sample_ID}_unmapped_1/2.fastq.gz)
        output_file (str): Path to write alignment summary TSV
        tmp_dir (str): Temp directory for intermediates
        run_alignment (bool): If False, skip alignments and return ""
        max_workers (int|None): Parallel jobs (default: half of logical cores, at least 1)
        bwa_threads_per_job (int): Threads per bwa mem process
        min_contig_len (int): Min contig length to align

    Returns:
        str: Path to TSV summary (or "" if run_alignment=False)
    """
    if not run_alignment:
        print("⚠️ Alignment disabled. Skipping.")
        return ""

    os.makedirs(tmp_dir, exist_ok=True)

    # Load contigs (filter by length)
    contig_dict: Dict[str, SeqIO.SeqRecord] = {
        rec.id: rec
        for rec in SeqIO.parse(merged_fasta, "fasta")
        if len(rec.seq) >= min_contig_len
    }
    print(f"🔍 Loaded {len(contig_dict)} contigs ≥ {min_contig_len} bp")

    # Load & sanitize hits
    df = pd.read_csv(diamond_tsv, sep="\t")
    # Normalize contig id column name
    if "Contig_ID" in df.columns:
        cid_col = "Contig_ID"
    elif "query_id" in df.columns:
        cid_col = "query_id"
    else:
        raise ValueError("Input TSV must contain 'Contig_ID' or 'query_id' column.")

    if "Sample_ID" not in df.columns:
        raise ValueError("Input TSV must contain 'Sample_ID' column.")

    # Keep only rows with contigs we have and de-duplicate (some DIAMOND rows may repeat the pair)
    pairs = (
        df[[ "Sample_ID", cid_col ]]
        .rename(columns={cid_col: "Contig_ID"})
        .dropna()
        .drop_duplicates()
    )

    # Pre-check FASTQs and filter
    def fq_paths(sample_id: str) -> Tuple[str, str]:
        r1 = os.path.join(fastq_dir, f"{sample_id}_unmapped_1.fastq.gz")
        r2 = os.path.join(fastq_dir, f"{sample_id}_unmapped_2.fastq.gz")
        return r1, r2

    mask_has_fastq = pairs["Sample_ID"].apply(lambda s: all(os.path.exists(p) for p in fq_paths(s)))
    missing_samples = pairs.loc[~mask_has_fastq, "Sample_ID"].unique().tolist()
    if missing_samples:
        print(f"[!] Missing FASTQs for {len(missing_samples)} sample(s): {missing_samples[:5]}{'...' if len(missing_samples)>5 else ''}")
    pairs = pairs[mask_has_fastq]

    # Filter to contigs we loaded
    mask_has_contig = pairs["Contig_ID"].isin(contig_dict.keys())
    missing_contigs = pairs.loc[~mask_has_contig, "Contig_ID"].unique().tolist()
    if missing_contigs:
        print(f"[!] {len(missing_contigs)} contig(s) not in merged_fasta or too short (skipping some jobs).")
    pairs = pairs[mask_has_contig].reset_index(drop=True)

    if pairs.empty:
        print("ℹ️ No valid (Sample_ID, Contig_ID) pairs to process after filtering.")
        # Still create an empty file with headers
        pd.DataFrame(columns=[
            "Sample_ID","Contig_ID","Contig_Length","Aligned_Reads","Breadth_Coverage","Mean_Depth"
        ]).to_csv(output_file, sep="\t", index=False)
        return output_file

    if max_workers is None:
        try:
            import os as _os
            logical = _os.cpu_count() or 2
        except Exception:
            logical = 2
        max_workers = max(1, logical // 2)

    print(f"🧵 Running up to {max_workers} jobs in parallel; BWA threads per job = {bwa_threads_per_job}")

    def process_one(sample_id: str, contig_id: str) -> Optional[Dict[str, Any]]:
        contig = contig_dict[contig_id]
        contig_len = len(contig.seq)

        # Construct per-job paths (safe filenames)
        safe_id = contig_id.replace("|", "_").replace("/", "_").replace("\\", "_").replace(" ", "_")
        prefix = f"{sample_id}_{safe_id}"
        ref_fasta = os.path.join(tmp_dir, f"{prefix}.fa")
        bam_raw = os.path.join(tmp_dir, f"{prefix}.bam")
        bam_sorted = os.path.join(tmp_dir, f"{prefix}.sorted.bam")

        r1, r2 = fq_paths(sample_id)

        # Write minimal single-contig reference
        try:
            if not os.path.exists(ref_fasta):
                SeqIO.write(contig, ref_fasta, "fasta")
        except Exception as e:
            print(f"[write fasta] {prefix}: {e}")
            return None

        # Index reference (once per job)
        try:
            needed = [".bwt", ".amb", ".ann", ".pac", ".sa"]
            if not all(os.path.exists(ref_fasta + ext) for ext in needed):
                subprocess.run(["bwa", "index", ref_fasta], check=True)
        except subprocess.CalledProcessError as e:
            print(f"[bwa index] {prefix}: {e}")
            _cleanup(prefix, ref_fasta, [bam_raw, bam_sorted])
            return None

        # Align if BAM absent
        if not os.path.exists(bam_sorted):
            try:
                if not os.path.exists(bam_raw):
                    bwa = subprocess.Popen(
                        ["bwa", "mem", "-t", str(bwa_threads_per_job), ref_fasta, r1, r2],
                        stdout=subprocess.PIPE
                    )
                    subprocess.run(
                        ["samtools", "view", "-bS", "-o", bam_raw],
                        stdin=bwa.stdout, check=True
                    )
                    bwa.stdout.close()
                    bwa.wait()

                # sort + index
                subprocess.run(["samtools", "sort", "-o", bam_sorted, bam_raw], check=True)
                subprocess.run(["samtools", "index", bam_sorted], check=True)
            except subprocess.CalledProcessError as e:
                print(f"[align/sort/index] {prefix}: {e}")
                _cleanup(prefix, ref_fasta, [bam_raw, bam_sorted])
                return None

        # Count mapped reads
        try:
            res = subprocess.run(
                ["samtools", "view", "-c", "-F", "4", bam_sorted],
                capture_output=True, text=True, check=True
            )
            mapped_reads = int(res.stdout.strip())
        except subprocess.CalledProcessError as e:
            print(f"[count mapped] {prefix}: {e}")
            mapped_reads = 0

        # Compute coverage & mean depth using samtools depth -aa
        # -aa ensures we include zero-coverage positions so mean depth uses contig_len as denominator.
        mean_depth = 0.0
        breadth_cov = 0.0
        try:
            # depth output: <contig>\t<pos>\t<depth>
            proc = subprocess.Popen(
                ["samtools", "depth", "-aa", bam_sorted],
                stdout=subprocess.PIPE, text=True
            )
            total_pos = 0
            covered_pos = 0
            sum_depth = 0

            # Because BAM has a single reference sequence, we can just iterate
            for line in proc.stdout:
                # Fast path: avoid split costs when line obviously too short
                if '\t' not in line:
                    continue
                parts = line.rstrip("\n").split("\t")
                # parts[0] is contig name in BAM (should match contig.id)
                # parts[1] is 1-based position
                # parts[2] is depth
                if len(parts) < 3:
                    continue
                try:
                    d = int(parts[2])
                except ValueError:
                    continue
                total_pos += 1
                if d > 0:
                    covered_pos += 1
                sum_depth += d

            proc.stdout.close()
            proc.wait()

            # Safety: if samtools didn't emit all positions, clamp denominators to contig_len
            denom = contig_len if contig_len > 0 else max(1, total_pos)
            if total_pos == 0:
                # Fallback to zero coverage if no output
                mean_depth = 0.0
                breadth_cov = 0.0
            else:
                # Use contig_len as denominator so zeros outside emitted rows are included
                mean_depth = float(sum_depth) / float(denom)
                breadth_cov = float(covered_pos) / float(denom)

        except Exception as e:
            print(f"[depth] {prefix}: {e}")
            mean_depth = 0.0
            breadth_cov = 0.0

        print(f"✅ {sample_id} vs {contig_id}: reads={mapped_reads}, cov={breadth_cov:.3f}, depth={mean_depth:.2f}")

        # Cleanup intermediates
        _cleanup(prefix, ref_fasta, [bam_raw, bam_sorted, bam_sorted + ".bai"])

        return {
            "Sample_ID": sample_id,
            "Contig_ID": contig_id,
            "Contig_Length": contig_len,
            "Aligned_Reads": mapped_reads,
            "Breadth_Coverage": round(breadth_cov, 6),
            "Mean_Depth": round(mean_depth, 6),
        }

    def _cleanup(prefix: str, ref_fasta: str, paths: list):
        # remove bam(s), bai, ref fasta & indexes
        for p in paths:
            try:
                if p and os.path.exists(p):
                    os.remove(p)
            except Exception:
                pass
        try:
            if os.path.exists(ref_fasta):
                os.remove(ref_fasta)
        except Exception:
            pass
        for ext in [".bwt", ".amb", ".ann", ".pac", ".sa"]:
            try:
                f = ref_fasta + ext
                if os.path.exists(f):
                    os.remove(f)
            except Exception:
                pass

    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        fut2key = {}
        for _, r in pairs.iterrows():
            fut = ex.submit(process_one, r["Sample_ID"], r["Contig_ID"])
            fut2key[fut] = (r["Sample_ID"], r["Contig_ID"])

        for fut in as_completed(fut2key):
            try:
                res = fut.result()
                if res:
                    results.append(res)
            except Exception as e:
                s, c = fut2key[fut]
                print(f"[job error] {s} vs {c}: {e}")

    # Save summary
    out_df = pd.DataFrame(results, columns=[
        "Sample_ID","Contig_ID","Contig_Length","Aligned_Reads","Breadth_Coverage","Mean_Depth"
    ])
    out_df.to_csv(output_file, sep="\t", index=False)
    print(f"\n📝 Alignment summary saved to: {output_file}")
    return output_file



