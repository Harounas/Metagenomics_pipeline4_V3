def save_genomad_short_virus_contigs(
    input_tsv: str,
    output_tsv: str = "short_contig_virus_table.tsv",
    max_virus_len: int = 500,   # virus kept only if its longest contig < 500 bp
    min_contig_len: int = 200,  # discard any contig shorter than this
    suffix_re: str = r"(.+?)_\d+$",  # strip trailing “_1”, “_2”, …
):
    """
    From a geNomad merged-genes TSV create a 4-column table:
        Sample_ID   Sample_ID|gene   contig_length   taxname

    • A virus (sample × taxname pair) is included only if *all* its contigs
      are < `max_virus_len` bp.
    • Individual contigs shorter than `min_contig_len` are discarded.
    • Trailing “_N” suffixes are removed from the contig part of the gene.
    • Underscores in taxname are replaced with spaces.
    """

    suffix_pat      = re.compile(suffix_re)
    contig_len_pat  = re.compile(r"length_(\d+)")      # pulls “230” from “…length_230_cov_…”
    # 1️⃣  First pass – collect max-contig length per (sample, taxname)
    max_len = defaultdict(int)   # (sample_id, taxname) ➜ longest contig length

    rows_tmp = []                # will hold rows for second pass
    with open(input_tsv, newline="") as fin:
        reader = csv.DictReader(fin, delimiter="\t")
        for row in reader:
            raw_gene  = row["gene"]
            clean_gene = suffix_pat.sub(r"\1", raw_gene)

            sample_id  = clean_gene.split("|", 1)[0]
            match_len  = contig_len_pat.search(clean_gene)
            if not match_len:
                continue                             # skip if length not found in ID
            contig_len = int(match_len.group(1))

            taxname = row.get("taxname", "").replace("_", " ")
            key = (sample_id, taxname)
            if contig_len > max_len[key]:
                max_len[key] = contig_len

            rows_tmp.append((sample_id, clean_gene, contig_len, taxname))

    # 2️⃣  Second pass – write rows that satisfy both filters
    with open(output_tsv, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(["Sample_ID", "Sample_ID|gene", "contig_length", "taxname"])

        kept = 0
        for sample_id, clean_gene, contig_len, taxname in rows_tmp:
            if contig_len < min_contig_len:
                continue                                            # too short
            if max_len[(sample_id, taxname)] >= max_virus_len:
                continue                                            # virus too long
            writer.writerow([sample_id, clean_gene, contig_len, taxname])
            kept += 1

    print(f"✓ Wrote {kept:,} contig rows to {output_tsv}")





# -----------------------------------------------------------
# 1) geNomad gene table → short-virus contigs (4-column TSV)
# -----------------------------------------------------------

def save_genomad_short_virus_contigs(
    input_tsv: str,
    output_tsv: str = "short_contig_virus_table.tsv",
    output_dir: str = ".",                 # NEW
    max_virus_len: int = 500,              # virus kept only if its longest contig < 500 bp
    min_contig_len: int = 200,             # discard any contig shorter than this
    suffix_re: str = r"(.+?)_\d+$"         # strip trailing “_1”, “_2”, …
):
    """
    Writes <output_dir>/<output_tsv> with columns:
        Sample_ID   Sample_ID|gene   contig_length   taxname
    Filters:
      • keep a virus (sample × taxname) only if *all* its contigs < max_virus_len
      • keep individual contigs only if length ≥ min_contig_len
    """
    suffix_pat      = re.compile(suffix_re)
    contig_len_pat  = re.compile(r"length_(\d+)")

    # ── 1️⃣ first pass: find longest contig per (sample, taxname)
    max_len = defaultdict(int)          # (sample, taxname) → longest length
    rows_tmp = []                       # cache rows for second pass

    with open(input_tsv, newline="") as fin:
        reader = csv.DictReader(fin, delimiter="\t")
        for row in reader:
            raw_gene   = row["gene"]
            clean_gene = suffix_pat.sub(r"\1", raw_gene)
            sample_id  = clean_gene.split("|", 1)[0]

            m_len = contig_len_pat.search(clean_gene)
            if not m_len:
                continue
            contig_len = int(m_len.group(1))

            taxname = row.get("taxname", "").replace("_", " ")
            key = (sample_id, taxname)
            max_len[key] = max(max_len[key], contig_len)

            rows_tmp.append((sample_id, clean_gene, contig_len, taxname))

    # ── 2️⃣ second pass: write filtered rows
    output_path = Path(output_dir) / output_tsv
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", newline="") as fout:
        w = csv.writer(fout, delimiter="\t")
        w.writerow(["Sample_ID", "Sample_ID|gene", "contig_length", "taxname"])

        kept = 0
        for sample_id, clean_gene, contig_len, taxname in rows_tmp:
            if contig_len < min_contig_len:
                continue
            if max_len[(sample_id, taxname)] >= max_virus_len:
                continue
            w.writerow([sample_id, clean_gene, contig_len, taxname])
            kept += 1

    print(f"✓ Wrote {kept:,} contig rows to {output_path}")



# -----------------------------------------------------------
# 2) Kraken-based filter of short-virus contigs
# -----------------------------------------------------------

def extract_short_contigs_kraken(
    base_contigs_dir,
    output_tsv="short_contigs_summary.tsv",
    output_dir: str = ".",                # NEW
    max_virus_len=500,                    # virus kept only if its longest contig < 500 bp
    min_contig_len=200,                   # discard any contig shorter than this
):
    """
    Writes <output_dir>/<output_tsv> with columns:
        Sample_ID   Sample_ID|gene   contig_length   taxname
    """
    allowed_ranks = {"S", "S1", "S2"}
    base_dir = Path(base_contigs_dir)
    output_path = Path(output_dir) / output_tsv
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", newline="") as out_f:
        writer = csv.writer(out_f, delimiter="\t")
        writer.writerow(["Sample_ID", "Sample_ID|gene", "contig_length", "taxname"])

        for report in base_dir.glob("*_Viruses_kraken_report.txt"):
            sample_id  = report.stem.replace("_Viruses_kraken_report", "")
            kout       = base_dir / f"{sample_id}_kraken2_output.txt"
            contigs_fa = base_dir / sample_id / "contigs.fasta"
            if not (kout.exists() and contigs_fa.exists()):
                print(f"⚠  Skipping {sample_id}: missing files.")
                continue

            # 1️⃣ taxid → virus name
            taxon_map = {}
            with open(report) as rf:
                for row in csv.reader(rf, delimiter="\t"):
                    if len(row) < 6 or row[3] not in allowed_ranks:
                        continue
                    taxid = row[4].strip()
                    name  = row[5]
                    if any(tok in name.lower() for tok in ("virus", "virinae", "viridae")):
                        taxon_map[taxid] = name
            if not taxon_map:
                continue

            # 2️⃣ taxid → contig IDs
            contig_hits = {}
            with open(kout) as kf:
                for row in csv.reader(kf, delimiter="\t"):
                    if len(row) < 3:
                        continue
                    cid, info = row[1].strip(), row[2]
                    for tid in taxon_map:
                        if f"taxid {tid}" in info:
                            contig_hits.setdefault(tid, []).append(cid)
            if not contig_hits:
                continue

            # 3️⃣ contig length lookup
            length_map = {rec.id: len(rec.seq) for rec in SeqIO.parse(contigs_fa, "fasta")}

            # 4️⃣ filter & write
            for tid, cids in contig_hits.items():
                lengths = [length_map.get(cid, 0) for cid in cids]
                if not lengths or max(lengths) >= max_virus_len:
                    continue

                vname = taxon_map[tid]
                for cid, clen in zip(cids, lengths):
                    if clen >= min_contig_len:
                        writer.writerow([sample_id, f"{sample_id}|{cid}", clen, vname])

    print(f"✓ Summary written to {output_path}")


def merge_short_contigs(
    tsv_a: str,
    tsv_b: str,
    output_tsv: str = "combined_no_dupes.tsv",
    output_dir: str = "."
):
    """
    Combine two TSV tables and remove duplicate (Sample_ID|gene, taxname) pairs.

    Parameters:
    -----------
    tsv_a : str
        Path to the first TSV file.
    tsv_b : str
        Path to the second TSV file.
    output_tsv : str
        Name of the output TSV file.
    output_dir : str
        Directory where output TSV will be written.

    Returns:
    --------
    Writes the combined TSV file to: {output_dir}/{output_tsv}
    """
    output_path = Path(output_dir) / output_tsv
    output_path.parent.mkdir(parents=True, exist_ok=True)  # ensure directory exists

    uniq = set()  # Track (Sample_ID|gene, taxname) pairs
    written = 0

    with open(output_path, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(["Sample_ID", "Sample_ID|gene", "contig_length", "taxname"])

        for src in (tsv_a, tsv_b):
            with open(src, newline="") as fin:
                reader = csv.reader(fin, delimiter="\t")
                header = next(reader)  # skip header
                for row in reader:
                    contig_id = row[1]
                    taxname   = row[3]
                    key = (contig_id, taxname)
                    if key in uniq:
                        continue
                    uniq.add(key)
                    writer.writerow(row)
                    written += 1

    print(f"✓ Merged:\n   {Path(tsv_a)}\n   {Path(tsv_b)}")
    print(f"→ Wrote {written:,} unique rows to {output_path}")


def run_diamond(diamond_db: str,
                query_file: str,
                output_file: str = "results.m8",
                threads: int = 8) -> None:
    """
    Runs Diamond BLASTX of query_file against diamond_db.
    """
    cmd = [
        "diamond", "blastx",
        "--query", query_file,
        "--db", diamond_db,
        "--out", output_file,
        "--threads", str(threads),
        "--outfmt", "6"  # tabular
    ]
    print("\n🚀 Executing Diamond BLASTX:\n", " ".join(cmd))
    subprocess.run(cmd, check=True)
    print(f"✅ Diamond results written to {output_file}")


