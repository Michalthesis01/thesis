#!/usr/bin/env python3
# filter_fasta_by_mismatch.py
# Filter a FASTA file based on best-alignment mismatch rate from a PAF file.
"""
FASTA Filtering by Best-Hit Mismatch Rate
=========================================

Description:
    This script executes a final, strict quality filtering step on a FASTA file using a single 
    mismatch rate cutoff derived from PAF alignment data. It identifies the single best alignment 
    for every query sequence and retains only those sequences whose best hit's mismatch rate 
    is below a specified threshold.

Key Features:
    1. Best-Hit Selection: For sequences with multiple alignments, it uses a stringent ranking 
       hierarchy to determine the single "best" hit per query: 
       **Primary Status $\to$ Alignment Length $\to$ Mismatch Rate $\to$ MAPQ.**
    2. Mismatch Rate Calculation: Accurately calculates the mismatch rate for the best hit, 
       prioritizing PAF tags (`de:f`, `NM:i`) over basic column arithmetic.
    3. Mismatch Filtering: Sequences are **failed** if they have no alignment in the PAF, or 
       if their best alignment's mismatch rate is strictly greater than the user-defined 
       `MISMATCH_THRESHOLD`.
    4. Output Categorization: Splits the input FASTA into two distinct output files: one 
       containing sequences that **passed** the mismatch filter, and one containing sequences 
       that **failed**.
    5. Quantitative Summary: Provides a detailed breakdown of the final pass/fail results, 
       quantifying both the count and total length (Mb and percentage) of retained/discarded sequences.
    6. Failing Query Report: Generates a TSV detailing the specific characteristics (length, 
       mismatch rate, MAPQ) of the best alignment for all queries that failed the threshold.

Parameters (User Inputs):
    PAF_PATH (str)             : Input PAF file containing all alignments.
    QUERY_FASTA (str)          : Input FASTA file (e.g., previously filtered CORE sequences).
    OUT_DIR (str)              : Directory for saving all outputs.
    MISMATCH_THRESHOLD (float) : The maximum allowable mismatch percentage for the best hit to pass.

Outputs:
    1. Filtered FASTA (PASSED): Sequences meeting the mismatch threshold.
    2. Filtered FASTA (FAILED): Sequences exceeding the threshold or unmapped.
    3. TSV: `failing_queries_...tsv` (Details of the best hit for all failed queries).
    4. TSV: `summary_pass_fail_lengths.tsv` (Overall Mb and percentage summary).
"""
import os, gzip
import pandas as pd
import numpy as np

# -------------------- USER PATHS --------------------
PAF_PATH = "/Users/michal/Desktop/pythonProject1_2/cut_genomes/alignments_with_cut_genomes/alignment_proper_cuts/filtered_alignments/filtered_fasta/hapil_to_florida_core_25kb_29_10.paf"
QUERY_FASTA = "/Users/michal/Desktop/pythonProject1_2/cut_genomes/alignments_with_cut_genomes/alignment_proper_cuts/filtered_alignments/filtered_fasta/filtered_fasta_25kb/hapil_25kb_CORE_only.fasta"
OUT_DIR = "/Users/michal/Desktop/pythonProject1_2/cut_genomes/alignments_with_cut_genomes/alignment_proper_cuts/filtered_alignments/filtered_fasta/quality_alignment_different_sorting_primary_length>mismatch>mapq/filtered_fasta"

MISMATCH_THRESHOLD = 60.0  # percent
# ----------------------------------------------------

os.makedirs(OUT_DIR, exist_ok=True)

def opengz(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")

def parse_paf_compute_mismatch(paf_path):
    """Load PAF and compute mismatch rate (%)."""
    rows = []
    with opengz(paf_path) as f:
        for ln in f:
            if not ln or ln[0] == "#":
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            try:
                qn, qlen, qs, qe = p[0], int(p[1]), int(p[2]), int(p[3])
                tn, tlen, ts, te = p[5], int(p[6]), int(p[7]), int(p[8])
                nmatch, blen, mapq = int(p[9]), int(p[10]), int(p[11])
            except Exception:
                continue

            nm, de, tp = None, None, "?"
            for opt in p[12:]:
                if opt.startswith("NM:i:"):
                    try: nm = int(opt.split(":")[-1])
                    except: pass
                elif opt.startswith("de:f:"):
                    try: de = float(opt.split(":")[-1])
                    except: pass
                elif opt.startswith("tp:A:"):
                    tp = opt.split(":")[-1]

            if blen <= 0:
                mismatch = np.nan
            elif de is not None:
                mismatch = (1 - de) * 100.0
            elif nm is not None:
                mismatch = (nm / blen) * 100.0
            else:
                mismatch = ((blen - nmatch) / blen) * 100.0

            rows.append({
                "query": qn,
                "qlen": qlen,
                "aln_span": qe - qs,
                "mapq": mapq,
                "tp": tp,
                "mismatch_rate_pct": mismatch
            })
    return pd.DataFrame(rows)

def fasta_lengths(path):
    """Return dict {id: length}."""
    lens = {}
    with opengz(path) as f:
        name, length = None, 0
        for ln in f:
            if ln.startswith(">"):
                if name is not None:
                    lens[name] = length
                name = ln[1:].strip().split()[0]
                length = 0
            else:
                length += len(ln.strip())
        if name is not None:
            lens[name] = length
    return lens

def fasta_iter_ids_only(path):
    """Return list of sequence IDs (first word after '>')."""
    ids = []
    with opengz(path) as f:
        for ln in f:
            if ln.startswith(">"):
                ids.append(ln[1:].strip().split()[0])
    return ids

def write_fasta_subset(in_fa, out_fa, keep_ids):
    """Write subset of sequences to new FASTA."""
    n_kept = 0
    with opengz(in_fa) as fin, open(out_fa, "w") as fout:
        write = False
        for ln in fin:
            if ln.startswith(">"):
                qid = ln[1:].strip().split()[0]
                write = qid in keep_ids
                if write:
                    fout.write(ln)
                    n_kept += 1
            else:
                if write:
                    fout.write(ln)
    return n_kept

# --- MAIN ---
print(f"[INFO] Loading PAF: {PAF_PATH}")
df = parse_paf_compute_mismatch(PAF_PATH)
if df.empty:
    raise SystemExit("[ERROR] No alignments found.")

# Sorting priority: primary first → longest → lowest mismatch → highest MAPQ
tp_priority = {"P": 0, "S": 1, "?": 2}
df["tp_rank"] = df["tp"].map(tp_priority).fillna(2)
df_sorted = df.sort_values(
    ["query", "tp_rank", "aln_span", "mismatch_rate_pct", "mapq"],
    ascending=[True, True, False, True, False]
)
best = df_sorted.groupby("query", as_index=False).first()

# FASTA lengths
print("[INFO] Reading FASTA lengths...")
fasta_lens = fasta_lengths(QUERY_FASTA)
fasta_ids = list(fasta_lens.keys())
fasta_id_set = set(fasta_ids)

# Build pass/fail sets
best_by_query = best.set_index("query", drop=False)
pass_ids, fail_ids = set(), set()

for q in fasta_ids:
    if q not in best_by_query.index:
        fail_ids.add(q)
        continue
    val = best_by_query.at[q, "mismatch_rate_pct"]
    if pd.isna(val) or float(val) > MISMATCH_THRESHOLD:
        fail_ids.add(q)
    else:
        pass_ids.add(q)

# Write outputs
fail_fa = os.path.join(
    OUT_DIR, f"{os.path.splitext(os.path.basename(QUERY_FASTA))[0]}__FAILED_mismatch_gt_{int(MISMATCH_THRESHOLD)}pct_or_unmapped.fasta"
)
pass_fa = os.path.join(
    OUT_DIR, f"{os.path.splitext(os.path.basename(QUERY_FASTA))[0]}__PASSED_mismatch_le_{int(MISMATCH_THRESHOLD)}pct.fasta"
)

n_fail = write_fasta_subset(QUERY_FASTA, fail_fa, fail_ids)
n_pass = write_fasta_subset(QUERY_FASTA, pass_fa, pass_ids)

# Compute total lengths (Mb)
len_pass_bp = sum(fasta_lens[q] for q in pass_ids if q in fasta_lens)
len_fail_bp = sum(fasta_lens[q] for q in fail_ids if q in fasta_lens)
len_total_bp = len_pass_bp + len_fail_bp
len_pass_mb = len_pass_bp / 1e6
len_fail_mb = len_fail_bp / 1e6
len_total_mb = len_total_bp / 1e6
pct_pass = (len_pass_bp / len_total_bp * 100) if len_total_bp else 0
pct_fail = (len_fail_bp / len_total_bp * 100) if len_total_bp else 0

# Write TSV of failing queries
fail_records = []
for q in sorted(fail_ids):
    if q in best_by_query.index:
        row = best_by_query.loc[q]
        fail_records.append({
            "query": q,
            "qlen": row.get("qlen", ""),
            "aln_span": row.get("aln_span", ""),
            "mismatch_rate_pct": row.get("mismatch_rate_pct", ""),
            "mapq": row.get("mapq", ""),
            "tp": row.get("tp", "")
        })
    else:
        fail_records.append({"query": q, "qlen": "", "aln_span": "", "mismatch_rate_pct": "", "mapq": "", "tp": ""})

fail_tsv = os.path.join(
    OUT_DIR, f"failing_queries_mismatch_gt_{int(MISMATCH_THRESHOLD)}pct_or_unmapped.tsv"
)
pd.DataFrame(fail_records).to_csv(fail_tsv, sep="\t", index=False)

# Write length summary TSV
summary_tsv = os.path.join(OUT_DIR, f"summary_pass_fail_lengths.tsv")
with open(summary_tsv, "w") as f:
    f.write("category\tcount\tlength_Mb\tpercent_of_total\n")
    f.write(f"PASSED\t{len(pass_ids)}\t{len_pass_mb:.3f}\t{pct_pass:.2f}\n")
    f.write(f"FAILED\t{len(fail_ids)}\t{len_fail_mb:.3f}\t{pct_fail:.2f}\n")
    f.write(f"TOTAL\t{len(fasta_ids)}\t{len_total_mb:.3f}\t100.00\n")

# Summary
print("\n=== FILTER SUMMARY ===")
print(f"Total queries in FASTA: {len(fasta_ids):,}")
print(f"Aligned (appear in PAF): {len(set(best['query'])):,}")
print(f"PASSED (≤ {MISMATCH_THRESHOLD}% mismatch): {len(pass_ids):,}  ({len_pass_mb:.2f} Mb, {pct_pass:.2f}%)  → {pass_fa}")
print(f"FAILED  (> {MISMATCH_THRESHOLD}% or unmapped): {len(fail_ids):,}  ({len_fail_mb:.2f} Mb, {pct_fail:.2f}%)  → {fail_fa}")
print(f"TOTAL genome length: {len_total_mb:.2f} Mb")
print(f"Summary written to: {summary_tsv}")
print("========================\n")