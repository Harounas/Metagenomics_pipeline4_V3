import os
import subprocess

def run_fastp(forward, reverse, base_name, output_dir, threads):
    """
    Runs fastp for quality trimming of FASTQ files.

    Returns:
        Tuple: (trimmed_forward, trimmed_reverse) if paired-end, else (trimmed_forward, None)
    """
    trimmed_forward = os.path.join(output_dir, f"{base_name}_trimmed_R1.fastq.gz")
    trimmed_reverse = os.path.join(output_dir, f"{base_name}_trimmed_R2.fastq.gz")
    html_report = os.path.join(output_dir, f"{base_name}_fastp.html")
    json_report = os.path.join(output_dir, f"{base_name}_fastp.json")

    fastp_cmd = [
        "fastp",
        "-i", forward,
        "-o", trimmed_forward,
        "-w", str(threads),
        "--html", html_report,
        "--json", json_report
    ]

    if reverse:
        fastp_cmd.extend([
            "-I", reverse,
            "-O", trimmed_reverse
        ])

    logging.info("Running fastp: " + " ".join(fastp_cmd))
    subprocess.run(fastp_cmd, check=True)

    return trimmed_forward, trimmed_reverse if reverse else None
