#!/usr/bin/env python3
"""
PAF Alignment Multi-Configuration Filtering Pipeline
====================================================

Description:
    This is the primary script for performing complex, multi-criteria filtering on Minimap2 
    PAF alignment files (typically from chunked genomes). It systematically applies a sweep 
    of per-alignment coverage and top-N alignment limits to generate detailed summaries of 
    which query sequences (chunks) pass or fail the final quality gates.

Key Features:
    1. Multi-Dimensional Sweep: Iterates through defined combinations of **minimum per-alignment 
       coverage** (e.g., 20%, 30%, 50%) and **top-N alignment mode** ('top5' vs. 'all').
    2. Alignment Selection: For each query, alignments are first filtered by the minimum 
       per-alignment coverage (`query_coverage > per_align_cov_min`). The remaining alignments 
       are then restricted based on the `top_n` mode.
    3. Multi-Criteria Pass/Fail: A query is marked 'passed' only if **both** conditions are met:
        a. **Mismatch Gate:** At least one selected alignment has a mismatch rate $\le$ `MISMATCH_RATE_MAX`.
        b. **Coverage Gate:** The **merged, unique sequence coverage** of all passing alignments 
           meets or exceeds the `TOTAL_COVERAGE_PASS` threshold.
    4. Interval Merging: Uses a helper function (`merge_intervals`) to accurately calculate the 
       unique base pairs covered by the passing alignments.
    5. Output Generation: Creates dedicated, structured output directories for each configuration, 
       containing two primary files per input PAF:
        * **Per-Query CSV:** Detailed record for every query sequence, including merged coverage, 
          thresholds, and the final 'passed_both' status.
        * **Summary TXT:** A concise summary of total Megabases (Mb) and percentage of the 
          query genome retained and discarded, used directly by subsequent plotting scripts.

Parameters (User Inputs):
    BASE_ALIGN_DIR (str)       : Directory containing the raw input .paf files.
    BASE_OUT_DIR (str)         : Root directory for creating all configuration subfolders.
    MISMATCH_RATE_MAX (float)  : Maximum allowed mismatch percentage for an alignment to be considered.
    TOTAL_COVERAGE_PASS (float): Minimum unique query coverage percentage required to pass the final gate.
    MIN_COVERAGES (list)       : Per-alignment coverage thresholds to iterate over.
    TOPN_MODES (list)          : '5' (Top 5 alignments) or 'all' (All alignments) to iterate over.

Outputs:
    - Multiple Directories (e.g., `filtered_alignments_min.coverage=50_mismatch=30%_topN=5/`).
    - Per-Query CSVs (e.g., `aln_...__coverage_summary_mismatch30.csv`).
    - Summary TXT files (`..._summary.txt`) for plotting scripts.
"""
import os, glob
import pandas as pd
import numpy as np

# ---------- Paths that match your plotting script ----------
BASE_ALIGN_DIR = '/Users/michal/Desktop/pythonProject1/cut_genomes/alignments_with_cut_genomes/alignment_proper_cuts'
BASE_OUT_DIR   = '/Users/michal/Desktop/pythonProject1/cut_genomes/alignments_with_cut_genomes/alignment_proper_cuts/filtered_alignments'

# Fixed thresholds (can tweak)
MISMATCH_RATE_MAX = 30.0
TOTAL_COVERAGE_PASS = 80.0

MIN_COVERAGES = [50.0, 30.0, 20.0]
TOPN_MODES    = ['5', 'all']   # strings to match your folder names

# ---------- your helper funcs copied (trimmed) ----------
def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: x[0])
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start <= merged[-1][1] :
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged

def process_query_group(group, top_n, mismatch_rate_max, total_coverage_pass):
    q_length = group['query_length'].iloc[0]
    if top_n == 'all':
        selected = group
    else:
        selected = group.nlargest(int(top_n), 'query_coverage') if len(group) > 0 else group
    selected_passing = selected[selected['mismatch_rate'] <= mismatch_rate_max].copy()
    intervals = list(zip(selected_passing['query_start'], selected_passing['query_end']))
    merged_intervals = merge_intervals(intervals) if intervals else []
    if merged_intervals:
        unique_bp = sum(e - s for s, e in merged_intervals)
    else:
        unique_bp = 0
    total_query_coverage = (unique_bp / q_length) * 100 if q_length > 0 else 0.0
    coverage_pass = (total_query_coverage >= total_coverage_pass)
    mismatch_any_pass = (len(selected_passing) > 0)
    passed_both = 'passed' if (coverage_pass and mismatch_any_pass) else "didn't pass"
    return pd.Series({
        'merged_query_start': min((s for s, e in merged_intervals), default=np.nan),
        'merged_query_end': max((e for s, e in merged_intervals), default=np.nan),
        'query_length': q_length,
        'total_query_coverage': round(total_query_coverage, 2),
        'coverage_threshold': total_coverage_pass,
        'coverage_pass': 'passed' if coverage_pass else "didn't pass",
        'mismatch_rate_threshold': mismatch_rate_max,
        'alignments_used_for_merging': len(selected_passing),
        'top_n_considered': len(selected),
        'passed_both': passed_both
    })

def process_one_paf(paf_path, out_csv_path, per_align_cov_min, top_n):
    # Load all alignments
    df = pd.read_csv(paf_path, sep='\t', header=None, usecols=range(12), engine='python')
    df.columns = [
        'query_name','query_length','query_start','query_end','strand',
        'target_name','target_length','target_start','target_end',
        'num_matches','alignment_length','mapq'
    ]

    # --- keep full query list for a fixed denominator ---
    all_queries = df[['query_name','query_length']].drop_duplicates().copy()

    # Per-alignment metrics
    df['query_coverage'] = ((df['query_end'] - df['query_start']) / df['query_length']) * 100
    df['mismatches']     = df['alignment_length'] - df['num_matches']
    df['mismatch_rate']  = np.where(
        df['alignment_length'] > 0,
        (df['mismatches'] / df['alignment_length']) * 100,
        np.nan
    )

    # Gate by per-alignment coverage
    df_f = df[df['query_coverage'] > per_align_cov_min].copy()

    if df_f.empty:
        # nobody has any passing alignments: output all queries as didn't pass
        out = all_queries.copy()
        out['merged_query_start'] = np.nan
        out['merged_query_end']   = np.nan
        out['total_query_coverage'] = 0.0
        out['coverage_threshold']   = TOTAL_COVERAGE_PASS
        out['coverage_pass']        = "didn't pass"
        out['mismatch_rate_threshold']   = MISMATCH_RATE_MAX
        out['alignments_used_for_merging'] = 0
        out['top_n_considered'] = 0
        out['passed_both'] = "didn't pass"
    else:
        # Per-query computation on filtered alignments
        out_pass = (
            df_f.groupby('query_name', as_index=False)
                .apply(lambda g: process_query_group(
                    g, top_n, MISMATCH_RATE_MAX, TOTAL_COVERAGE_PASS
                ), include_groups=False)
                .reset_index(drop=True)
        )

        # --- left-join onto all queries to keep the missing ones as "didn't pass" ---
        out = all_queries.merge(out_pass, on=['query_name','query_length'], how='left')

        # fill queries with no rows in df_f
        out['merged_query_start'] = out['merged_query_start'].astype('float64')
        out['merged_query_end']   = out['merged_query_end'].astype('float64')
        out['total_query_coverage'] = out['total_query_coverage'].fillna(0.0)
        out['coverage_threshold']   = out['coverage_threshold'].fillna(TOTAL_COVERAGE_PASS)
        out['coverage_pass']        = out['coverage_pass'].fillna("didn't pass")
        out['mismatch_rate_threshold'] = out['mismatch_rate_threshold'].fillna(MISMATCH_RATE_MAX)
        out['alignments_used_for_merging'] = out['alignments_used_for_merging'].fillna(0).astype(int)
        out['top_n_considered'] = out['top_n_considered'].fillna(0).astype(int)
        out['passed_both'] = out['passed_both'].fillna("didn't pass")

    # neatness
    out['total_query_coverage'] = out['total_query_coverage'].round(2)

    cols = [
        'query_name',
        'merged_query_start','merged_query_end','query_length',
        'total_query_coverage','coverage_threshold','coverage_pass',
        'mismatch_rate_threshold','alignments_used_for_merging','top_n_considered',
        'passed_both'
    ]
    out[cols].to_csv(out_csv_path, index=False, float_format='%.2f')

    # Return totals for summary writing (fixed denominator = all_queries)
    total_mb  = all_queries['query_length'].sum() / 1_000_000.0
    passed_mb = out.loc[out['passed_both']=='passed', 'query_length'].sum() / 1_000_000.0
    return total_mb, passed_mb


def write_summary_txt(csv_path, txt_path, total_mb, passed_mb):
    df = pd.read_csv(csv_path)
    counts = df['passed_both'].value_counts() if 'passed_both' in df else pd.Series()
    lengths = df.groupby('passed_both')['query_length'].sum() if 'query_length' in df else pd.Series()
    lengths_mb = (lengths / 1_000_000).round(2) if not lengths.empty else pd.Series()

    discarded_mb = round(max(total_mb - passed_mb, 0), 2)
    retained_pct  = (passed_mb / total_mb * 100) if total_mb else 0.0
    discarded_pct = (discarded_mb / total_mb * 100) if total_mb else 0.0

    lines = []
    lines.append(f"=== Summary of '{os.path.basename(csv_path)}' ===\n")
    for status in ['passed', "didn't pass"]:
        n = int(counts.get(status, 0))
        mb = float(lengths_mb.get(status, passed_mb if status=='passed' else discarded_mb) or 0.0)
        lines.append(f"{status:12}: {n:5d} queries, {mb:8.2f} Mb total")
    lines.append(f"\nTotal genome length represented: {total_mb:.2f} Mb")
    lines.append(f"\nGenome retained:  {retained_pct:.2f}%")
    lines.append(f"Genome discarded: {discarded_pct:.2f}%")

    with open(txt_path, 'w') as f:
        f.write("\n".join(lines))

def main():
    paf_paths = sorted(glob.glob(os.path.join(BASE_ALIGN_DIR, "*.paf")))
    if not paf_paths:
        print(f"[ERROR] No .paf in {BASE_ALIGN_DIR}")
        return

    for cov in MIN_COVERAGES:
        for topn in TOPN_MODES:
            cfg_name = f"filtered_alignments_min.coverage={int(cov)}_mismatch=30%_topN={topn}"
            out_dir = os.path.join(BASE_OUT_DIR, cfg_name)
            os.makedirs(out_dir, exist_ok=True)
            results_dir = os.path.join(out_dir, "results")
            os.makedirs(results_dir, exist_ok=True)

            print(f"\n=== Running config: min.cov={cov}, topN={topn} => {cfg_name}")
            for paf in paf_paths:
                stem = os.path.splitext(os.path.basename(paf))[0]
                csv_name = f"{stem}__coverage_summary_mismatch{int(MISMATCH_RATE_MAX)}.csv"
                csv_path = os.path.join(out_dir, csv_name)

                total_mb, passed_mb = process_one_paf(
                    paf_path=paf,
                    out_csv_path=csv_path,
                    per_align_cov_min=cov,
                    top_n=topn
                )

                # summary .txt next to CSVs but under results/ per your plotter
                txt_name = f"{os.path.splitext(os.path.basename(csv_path))[0]}_summary.txt"
                txt_path = os.path.join(results_dir, txt_name)
                write_summary_txt(csv_path, txt_path, total_mb, passed_mb)

            print(f"[DONE] Wrote CSVs to: {out_dir}\n       Summaries to: {results_dir}")

if __name__ == "__main__":
    main()
