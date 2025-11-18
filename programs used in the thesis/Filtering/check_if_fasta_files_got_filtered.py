#!/usr/bin/env python3
"""
FASTA Filtering Audit and Quantification
========================================

Description:
    This script performs a critical audit of the FASTA filtering process. It compares the 
    list of query sequences present in the newly generated "CORE" FASTA files against the 
    pass/fail status recorded in the original QC summary CSVs. Its primary goal is to 
    verify that the correct sequences were retained and to quantify the total Megabases (Mb) 
    and percentage of the original query genome that was successfully kept versus removed.

Key Features:
    1. Integrity Checks (Audit): Performs sanity checks to identify discrepancies:
        * Sequences present in the CORE FASTA but not marked 'passed' in the CSV.
        * Sequences marked 'passed' in the CSV but missing from the CORE FASTA.
    2. Sequence Length Quantification: Reads lengths from both the original and CORE FASTA 
       files, but uses the **original FASTA lengths** for the final Mb calculations to 
       maintain consistency.
    3. Final Metrics Calculation: Calculates and prints the total query genome length (Mb) 
       considered, the Mb and count of sequences **kept**, and the final **retention/removal 
       percentage** for each genome pair.

Parameters (User Inputs):
    FASTA_HAPIL_ORIG, FASTA_FLORIDA_ORIG (str): Paths to the initial, unfiltered FASTA files.
    FASTA_HAPIL_CORE, FASTA_FLORIDA_CORE (str): Paths to the filtered "CORE-only" FASTA files.
    CSV_HAPIL, CSV_FLORIDA (str): Paths to the CSV files containing the `passed_both` status.

Outputs:
    - Console Summary: Prints detailed statistics including total Mb, Mb kept/removed, 
      query counts, and final retention/removal percentages, along with any integrity warnings.
"""
import os
import sys
import pandas as pd

# -------------------- Inputs --------------------
# Original 25kb FASTAs (before filtering)
FASTA_HAPIL_ORIG = "/Users/michal/Desktop/pythonProject1/cut_genomes/hapil_20181217_consensus_25kb.fasta"
FASTA_FLORIDA_ORIG = "/Users/michal/Desktop/pythonProject1/cut_genomes/florida_brilliance_25kb.fasta"

# CORE-only FASTAs you just generated
FASTA_HAPIL_CORE = "/Users/michal/Desktop/pythonProject1/cut_genomes/alignments_with_cut_genomes/alignment_proper_cuts/filtered_alignments/filtered_fasta/hapil_20181217_consensus_25kb__CORE_only.fasta"
FASTA_FLORIDA_CORE = "/Users/michal/Desktop/pythonProject1/cut_genomes/alignments_with_cut_genomes/alignment_proper_cuts/filtered_alignments/filtered_fasta/florida_brilliance_25kb__CORE_only.fasta"

# CSVs with per-query pass/fail (topN=all, cov=20, mismatch=30, total=80)
CSV_HAPIL = "/Users/michal/Desktop/pythonProject1/cut_genomes/alignments_with_cut_genomes/alignment_proper_cuts/filtered_alignments/filtered_alignments_min.coverage=20_mismatch=30%_topN=all/aln_hapil_to_red_gauntlet_25kb__coverage_summary_mismatch30.csv"
CSV_FLORIDA = "/Users/michal/Desktop/pythonProject1/cut_genomes/alignments_with_cut_genomes/alignment_proper_cuts/filtered_alignments/filtered_alignments_min.coverage=20_mismatch=30%_topN=all/aln_florida_to_royal_royce_25kb__coverage_summary_mismatch30.csv"
# ------------------------------------------------

MB = 1_000_000.0

def load_passed_and_denominator(csv_path):
    """Return: passed_set, denom_set (all queries present in CSV, regardless of pass/fail)."""
    df = pd.read_csv(csv_path)
    if "query_name" not in df.columns or "passed_both" not in df.columns:
        raise ValueError(f"CSV missing required columns: {csv_path}")
    df["query_name"] = df["query_name"].astype(str)
    denom = set(df["query_name"])
    passed = set(df.loc[df["passed_both"].astype(str).str.lower() == "passed", "query_name"])
    return passed, denom

def fasta_lengths(fasta_path):
    """Stream a FASTA and return dict {header_first_token: sequence_length}."""
    lens = {}
    with open(fasta_path, "r") as f:
        name = None
        length = 0
        for line in f:
            if not line:
                continue
            if line.startswith(">"):
                # flush previous
                if name is not None:
                    lens[name] = lens.get(name, 0) + length
                header = line[1:].strip()
                name = header.split()[0]
                length = 0
            else:
                length += len(line.strip())
        # flush last
        if name is not None:
            lens[name] = lens.get(name, 0) + length
    return lens

def summarize_pair(label, orig_fa, core_fa, csv_path):
    # Load CSV sets
    passed_set, denom_set = load_passed_and_denominator(csv_path)

    # Load FASTA lengths
    orig_len = fasta_lengths(orig_fa)
    core_len = fasta_lengths(core_fa)

    # Denominator = queries present in CSV AND present in the original FASTA (for length)
    denom = denom_set.intersection(orig_len.keys())

    # Kept set = names present in CORE fasta ∩ denom
    kept = set(core_len.keys()).intersection(denom)
    removed = denom - kept

    # Sanity checks
    extra_core = set(core_len.keys()) - passed_set
    if extra_core:
        print(f"[WARN] {label}: {len(extra_core)} sequences in CORE FASTA not marked 'passed' in CSV (first few): {list(sorted(extra_core))[:5]}")

    missing_passed = passed_set.intersection(denom) - set(core_len.keys())
    if missing_passed:
        print(f"[WARN] {label}: {len(missing_passed)} 'passed' sequences (in CSV) not found in CORE FASTA (first few): {list(sorted(missing_passed))[:5]}")

    # Sum lengths using ORIGINAL FASTA lengths for consistency
    denom_bp  = sum(orig_len[q] for q in denom)
    kept_bp   = sum(orig_len[q] for q in kept)
    removed_bp = denom_bp - kept_bp

    denom_mb = denom_bp / MB
    kept_mb  = kept_bp / MB
    removed_mb = removed_bp / MB

    kept_pct = (kept_mb / denom_mb * 100.0) if denom_mb > 0 else 0.0
    removed_pct = 100.0 - kept_pct

    print(f"\n=== {label} (FASTA check) ===")
    print(f"Original FASTA: {os.path.basename(orig_fa)}   CORE FASTA: {os.path.basename(core_fa)}")
    print(f"Denominator (CSV∩origFASTA) queries: {len(denom):,}")
    print(f"Queries kept: {len(kept):,}   removed: {len(removed):,}")
    print(f"Genome (Mb) — total: {denom_mb:.2f}   kept: {kept_mb:.2f}   removed: {removed_mb:.2f}")
    print(f"Retention: {kept_pct:.2f}%   |   Removed: {removed_pct:.2f}%")

def main():
    # Hapil → Red Gauntlet
    for path in [FASTA_HAPIL_ORIG, FASTA_HAPIL_CORE, CSV_HAPIL]:
        if not os.path.isfile(path):
            sys.exit(f"Missing file: {path}")
    summarize_pair("Hapil→Red Gauntlet", FASTA_HAPIL_ORIG, FASTA_HAPIL_CORE, CSV_HAPIL)

    # Florida → Royal Royce
    for path in [FASTA_FLORIDA_ORIG, FASTA_FLORIDA_CORE, CSV_FLORIDA]:
        if not os.path.isfile(path):
            sys.exit(f"Missing file: {path}")
    summarize_pair("Florida Brilliance→Royal Royce", FASTA_FLORIDA_ORIG, FASTA_FLORIDA_CORE, CSV_FLORIDA)

if __name__ == "__main__":
    main()
