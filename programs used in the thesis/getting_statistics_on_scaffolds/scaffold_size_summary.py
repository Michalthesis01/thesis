#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Scaffold size summary & cut-size justification (MATCHES your splitter logic).

- Scans all .fasta/.fa (optionally .gz) in INPUT_DIR
- For each genome file, computes:
    * # scaffolds, total bp
    * mean, median, min, max scaffold length
    * N50/L50, N90/L90
    * GC% (A/T/G/C only), N% (ambiguous)
    * length-bin distribution: <1k, 1–25k, 25–50k, 50–100k, ≥100k (by bp share)
    * for CUT_SIZES = {1k, 25k, 50k, 100k} using your EXACT split rules:
        - scaffolds ≥ cut: count, % of scaffolds, % of total bp
        - estimated number of chunks after cutting (with tiny-tail merge rule)
        - total bp merged into tails (i.e., tails < min_tail_frac*cut that were absorbed)
        - total bp placed in short last-chunks (i.e., tails ≥ min_tail_frac*cut)
        - average resulting chunk length (post-cut, theoretical)
- Writes a single consolidated .txt report to OUTPUT_DIR/scaffold_size_summary.txt

Defaults :
  INPUT_DIR  = "/Users/michal/Desktop/thesis_2/data/genomes"
  OUTPUT_DIR = "/Users/michal/Desktop/pythonProject1_2/scaffold_size"
  CUT_SIZES  = [1_000, 25_000, 50_000, 100_000]
  MIN_TAIL_FRAC = 0.30  # tiny tails are merged (same as splitter)

Usage (optional overrides):
  python scaffold_size_summary.py --input <dir> --output <dir> --min-tail-frac 0.30
"""

import os
import gzip
import argparse
from glob import glob
from statistics import median

# ------- Defaults -------
DEFAULT_INPUT  = "/Users/michal/Desktop/thesis_2/data/genomes"
DEFAULT_OUTPUT = "/Users/michal/Desktop/pythonProject1_2/scaffold_size"

CUT_SIZES = [1_000, 25_000, 50_000, 100_000]  # 1kb, 25kb, 50kb, 100kb
MIN_TAIL_FRAC_DEFAULT = 0.30

BINS = [
    ("<1k",            0,      1_000),
    ("1–25k",          1_000,  25_000),
    ("25–50k",         25_000, 50_000),
    ("50–100k",        50_000, 100_000),
    ("≥100k",          100_000, float("inf")),
]

def open_maybe_gzip(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")

def parse_fasta_lengths_and_gc(fpath):
    lengths = []
    gc = at = n = 0
    cur_len = 0
    with open_maybe_gzip(fpath) as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if cur_len > 0:
                    lengths.append(cur_len)
                    cur_len = 0
                continue
            s = line.strip()
            if not s:
                continue
            cur_len += len(s)
            ls = s.lower()
            gc += ls.count("g") + ls.count("c")
            at += ls.count("a") + ls.count("t")
            n  += ls.count("n")
    if cur_len > 0:
        lengths.append(cur_len)
    return lengths, gc, at, n

def nxx_lengths(lengths, q=0.5):
    if not lengths:
        return 0, 0
    total = sum(lengths)
    target = total * q
    lens_sorted = sorted(lengths, reverse=True)
    cum = 0
    for i, L in enumerate(lens_sorted, start=1):
        cum += L
        if cum >= target:
            return L, i
    return lens_sorted[-1], len(lens_sorted)

def humanize_bp(x):
    if x >= 1_000_000_000:
        return f"{x/1_000_000_000:.2f} Gb"
    if x >= 1_000_000:
        return f"{x/1_000_000:.2f} Mb"
    if x >= 1_000:
        return f"{x/1_000:.2f} kb"
    return f"{x} bp"

def estimate_chunks_with_tail_rule(lengths, cut, min_tail_frac):
    """
    Implements EXACTLY your split_record_to_chunks():
      - full_chunks = L // cut
      - tail = L % cut
      - if tail == 0: total = full_chunks
        elif tail < min_tail_frac*cut and full_chunks > 0: total = full_chunks  (tail merged)
        else: total = full_chunks + 1  (last chunk is short but kept)
    Also returns:
      - merged_tails_bp: sum of tails that were merged
      - short_last_chunks_bp: sum of tails that formed their own last chunk
      - total_chunk_bp: theoretical total bp across all chunks (equals sum(lengths))
      - avg_chunk_len: total_bp / total_chunks
    """
    total_chunks = 0
    merged_tails_bp = 0
    short_last_chunks_bp = 0
    for L in lengths:
        if L == 0:
            continue
        full = L // cut
        tail = L % cut
        if tail == 0:
            total = full
        else:
            if tail < int(min_tail_frac * cut) and full > 0:
                total = full
                merged_tails_bp += tail
            else:
                total = full + 1
                short_last_chunks_bp += tail
        total_chunks += total
    total_bp = sum(lengths)
    avg_chunk_len = (total_bp / total_chunks) if total_chunks > 0 else 0.0
    return total_chunks, merged_tails_bp, short_last_chunks_bp, avg_chunk_len

def summarize_genome(fpath, min_tail_frac):
    lengths, gc, at, n = parse_fasta_lengths_and_gc(fpath)
    if not lengths:
        return {"file": os.path.basename(fpath), "empty": True}

    total_bp = sum(lengths)
    num = len(lengths)
    minL = min(lengths)
    maxL = max(lengths)
    meanL = total_bp / num
    medL = median(lengths)
    N50, L50 = nxx_lengths(lengths, 0.5)
    N90, L90 = nxx_lengths(lengths, 0.9)
    gc_percent = (gc / max(gc + at, 1)) * 100
    n_percent  = (n  / max(total_bp, 1)) * 100

    # Bins (by bp share)
    bin_stats = []
    for name, lo, hi in BINS:
        if hi == float("inf"):
            bp = sum(L for L in lengths if L >= lo)
        else:
            bp = sum(L for L in lengths if (L >= lo and L < hi))
        pct = (bp / total_bp) * 100 if total_bp else 0.0
        bin_stats.append((name, pct))

    # Cut-size table with your exact rules
    cut_rows = []
    for cut in CUT_SIZES:
        ge_cnt = sum(1 for L in lengths if L >= cut)
        ge_bp  = sum(L for L in lengths if L >= cut)
        ge_cnt_pct = (ge_cnt / num) * 100 if num else 0.0
        ge_bp_pct  = (ge_bp  / total_bp) * 100 if total_bp else 0.0

        est_chunks, merged_tails_bp, short_last_chunks_bp, avg_chunk_len = \
            estimate_chunks_with_tail_rule(lengths, cut, min_tail_frac)

        cut_rows.append({
            "cut": cut,
            "ge_cnt": ge_cnt,
            "ge_cnt_pct": ge_cnt_pct,
            "ge_bp_pct": ge_bp_pct,
            "est_chunks": est_chunks,
            "merged_tails_bp": merged_tails_bp,
            "short_last_chunks_bp": short_last_chunks_bp,
            "avg_chunk_len": avg_chunk_len
        })

    return {
        "file": os.path.basename(fpath),
        "empty": False,
        "num_scaffolds": num,
        "total_bp": total_bp,
        "min": minL,
        "max": maxL,
        "mean": meanL,
        "median": medL,
        "N50": N50,
        "L50": L50,
        "N90": N90,
        "L90": L90,
        "gc_percent": gc_percent,
        "n_percent": n_percent,
        "bins": bin_stats,
        "cuts": cut_rows,
    }

def write_report(results, out_dir, min_tail_frac):
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "scaffold_size_summary.txt")
    with open(out_path, "w") as w:
        w.write("# Scaffold Size Summary & Cut-Size Justification (exact splitter logic)\n")
        w.write(f"# min_tail_frac = {min_tail_frac:.2f}\n\n")
        for R in results:
            w.write(f"=== {R['file']} ===\n")
            if R.get("empty"):
                w.write("No sequences found.\n\n")
                continue
            w.write(f"Scaffolds:               {R['num_scaffolds']}\n")
            w.write(f"Total size:              {R['total_bp']} bp ({humanize_bp(R['total_bp'])})\n")
            w.write(f"Mean length:             {R['mean']:.2f} bp\n")
            w.write(f"Median length:           {R['median']:.2f} bp\n")
            w.write(f"Min length:              {R['min']} bp\n")
            w.write(f"Max length:              {R['max']} bp\n")
            w.write(f"N50 / L50:               {R['N50']} bp / {R['L50']} scaffolds\n")
            w.write(f"N90 / L90:               {R['N90']} bp / {R['L90']} scaffolds\n")
            w.write(f"GC% (A/T/G/C only):      {R['gc_percent']:.2f}%\n")
            w.write(f"N (ambiguous) content:   {R['n_percent']:.2f}% of total bp\n")
            w.write("\nLength bins (bp share):\n")
            for name, pct in R["bins"]:
                w.write(f"  {name:<9} {pct:>6.2f}%\n")
            w.write("\nCut-size justification (APPLIES your tail-merge rule):\n")
            w.write("  Cut   ≥Cut scaffolds     ≥Cut bp     Est. chunks   Merged tails bp   Short last-chunks bp   Avg chunk len\n")
            for row in R["cuts"]:
                w.write(
                    f"  {humanize_bp(row['cut']):<6} "
                    f"{row['ge_cnt']:>6} / {row['ge_cnt_pct']:>5.1f}%   "
                    f"{row['ge_bp_pct']:>6.1f}%   "
                    f"{row['est_chunks']:>11}   "
                    f"{row['merged_tails_bp']:>15}   "
                    f"{row['short_last_chunks_bp']:>20}   "
                    f"{row['avg_chunk_len']:>13.1f}\n"
                )
            w.write("\nNotes:\n")
            w.write("- 'Merged tails bp' sums tails (< min_tail_frac*cut) that were absorbed into the previous chunk.\n")
            w.write("- 'Short last-chunks bp' sums tails (≥ min_tail_frac*cut) that became their own final chunk.\n")
            w.write("- 'Avg chunk len' reflects the expected average fragment length after cutting at that size.\n")
            w.write("- These numbers directly justify your chosen window sizes by quantifying fragmentation vs. efficiency.\n\n\n")
    return out_path

def main():
    ap = argparse.ArgumentParser(description="Summarize scaffold sizes and cut-size justification (exact splitter logic).")
    ap.add_argument("--input",  default=DEFAULT_INPUT,  help="Directory with .fasta/.fa (.gz ok).")
    ap.add_argument("--output", default=DEFAULT_OUTPUT, help="Directory to write scaffold_size_summary.txt.")
    ap.add_argument("--min-tail-frac", type=float, default=MIN_TAIL_FRAC_DEFAULT,
                    help="Tiny tail threshold; tails < fraction*cut are merged (default 0.30).")
    args = ap.parse_args()

    patterns = ["*.fasta", "*.fa", "*.fasta.gz", "*.fa.gz"]
    files = []
    for pat in patterns:
        files.extend(glob(os.path.join(args.input, pat)))
    files.sort()

    os.makedirs(args.output, exist_ok=True)

    if not files:
        out_path = os.path.join(args.output, "scaffold_size_summary.txt")
        with open(out_path, "w") as w:
            w.write("# No FASTA files found.\n")
            w.write(f"# Checked: {args.input}\n")
        print(f"Wrote empty report: {out_path}")
        return

    results = [summarize_genome(f, args["min_tail_frac"] if isinstance(args, dict) else args.min_tail_frac) for f in files]
    out_path = write_report(results, args.output, args.min_tail_frac)
    print(f"Wrote report: {out_path}")

if __name__ == "__main__":
    main()
