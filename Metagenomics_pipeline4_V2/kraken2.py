import os
import subprocess

def run_kraken2(forward, reverse, base_name, kraken_db, output_dir, threads):
    kraken_report = os.path.join(output_dir, f"{base_name}_kraken_report.txt")
    kraken_output = os.path.join(output_dir, f"{base_name}_kraken2_output.txt")

    # Check if we are using contigs.fasta (i.e., reverse is None)
    if reverse is None:
        # Use Kraken2 with contigs file (single file input)
        kraken_cmd = [
            "kraken2", "--db", kraken_db,
            "--report", kraken_report,
            "--output", kraken_output,
            "--threads", str(threads),
            "--use-names",
            forward  # The contigs file (forward)
        ]
    else:
        # Standard paired-end run
        kraken_cmd = [
            "kraken2", "--db", kraken_db,
            "--report", kraken_report,
            "--output", kraken_output,
            "--threads", str(threads),
            "--use-names",
            "--paired", "--gzip-compressed", forward, reverse
        ]

    # Print the Kraken2 command for debugging
    print("Running Kraken2 command:", " ".join(kraken_cmd))

    # Execute the Kraken2 command
    subprocess.run(kraken_cmd, check=True)

    return kraken_report
