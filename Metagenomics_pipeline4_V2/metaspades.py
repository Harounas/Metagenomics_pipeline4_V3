import os
import subprocess
import logging

def run_spades(forward, reverse, base_name, output_dir, threads=8):
    """
    Runs MetaSPAdes for de novo assembly.

    Parameters:
    - forward (str): Path to forward reads.
    - reverse (str): Path to reverse reads.
    - base_name (str): Sample identifier.
    - output_dir (str): Directory for output files.
    - threads (int): Number of CPU threads.

    Returns:
    - str: Path to the assembled contigs file.
    """
    sample_outdir = os.path.join(output_dir, base_name)
    contigs_output = os.path.join(sample_outdir, "contigs.fasta")

    if os.path.exists(contigs_output):
        logging.info(f"[SKIP] SPAdes assembly already exists for {base_name}.")
        return contigs_output

    os.makedirs(sample_outdir, exist_ok=True)

    cmd = [
        "metaspades.py",
        "-1", forward,
        "-2", reverse,
        "-o", sample_outdir,
        "-t", str(threads)
    ]

    logging.info(f"[RUN] Running MetaSPAdes for {base_name}: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"SPAdes failed for {base_name} (exit {e.returncode})") from e

    if os.path.exists(contigs_output):
        logging.info(f"[DONE] Assembly complete for {base_name}: {contigs_output}")
        return contigs_output
    else:
        raise RuntimeError(f"SPAdes finished but no contigs.fasta found for {base_name}")

