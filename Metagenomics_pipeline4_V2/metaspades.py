import os

def run_spades(forward, reverse, base_name, output_dir, threads):
    """
    Runs SPAdes for de novo assembly.

    Parameters:
    - forward (str): Path to forward reads.
    - reverse (str): Path to reverse reads.
    - base_name (str): Sample identifier.
    - output_dir (str): Directory for output files.
    - threads (int): Number of CPU threads.

    Returns:
    - str: Path to the assembled contigs file.
    """
    contigs_output = os.path.join(output_dir, base_name, "contigs.fasta")
    if os.path.exists(contigs_output):
        print(f"SPAdes assembly already exists for {base_name}. Using existing contigs.")
        return contigs_output

    spades_cmd = f"metaspades.py -1 {forward} -2 {reverse} -o {os.path.join(output_dir, base_name)} -t {threads}"
    os.system(spades_cmd)

    if os.path.exists(contigs_output):
        print(f"Assembly complete for {base_name}. Contigs saved at {contigs_output}")
        return contigs_output
    else:
        raise RuntimeError(f"SPAdes failed for {base_name}. No contigs.fasta found.")
