#!/usr/bin/env python3
"""
Uncut Genome Alignment Retention Comparison
===========================================

Description:
    This script is a post-processing tool that analyzes summary log files (from a prior 
    filtering step) that pertain specifically to **uncut** (full-length) genome alignment 
    results. It extracts retention statistics, aggregates them by filtering configuration 
    (min coverage and top mode), and generates a comparative bar chart.

Key Features:
    1. Configuration Extraction: Uses Regex on filenames to extract key filtering parameters 
       like 'min_coverage' (e.g., 20, 50%) and 'top_mode' ('Top 5' vs. 'All alignments').
    2. Robust Metrics Parsing: Reads internal metrics (Mb retained/discarded, retention %) 
       from the summary text files. Includes logic to infer missing percentages or MB values 
       if partial data is present in the logs.
    3. Comparative Plotting: Generates a multi-panel bar chart (one panel per genome pair) 
       to visualize the **retained genome percentage** as a function of the minimum coverage 
       threshold.
    4. Top Mode Comparison: Within each plot, results are grouped to directly compare the 
       effect of keeping only the 'Top 5' alignments versus considering 'All alignments'.

Parameters (User Inputs):
    RESULTS_DIR (str): Path to the folder containing the input summary .txt files.
    OUT_DIR (str)    : Directory where the output CSV and plot will be saved.
    PAIR_LABELS (dict): Dictionary mapping filename prefixes to human-readable genome pair labels.

Outputs:
    1. CSV: 'uncut_retained_discarded_summary.csv' (Tabular data of all parsed metrics).
    2. Plot: 'uncut_retained_by_min_coverage_top5_vs_all.png' (Comparative bar chart).
"""
import os
import re
import glob
import math
import pandas as pd
import matplotlib.pyplot as plt

# -------------------- Paths --------------------
RESULTS_DIR = "/Users/michal/Desktop/pythonProject1_2/filtered_alignments/results"
OUT_DIR     = os.path.join(os.path.dirname(RESULTS_DIR), "plots_uncut")
os.makedirs(OUT_DIR, exist_ok=True)
# ------------------------------------------------

# Recognize pairs from filename prefixes
PAIR_LABELS = {
    'fixed_florida_to_royal_royce': 'Florida Brilliance→Royal Royce',
    'fixed_hapil_to_red_gauntlet': 'Hapil→Red Gauntlet',
}

# Filename pattern pieces:
# ..._top5_cov20_mismatch30_total80_summary.txt
# ..._topall_cov50_mismatch30_total80_summary.txt
RE_TOP   = re.compile(r'_top(5|all)_', re.IGNORECASE)
RE_COV   = re.compile(r'_cov(\d+)_', re.IGNORECASE)

# Robust parsers for the .txt contents
RE_RETAINED     = re.compile(r"Genome retained:\s*([0-9.]+)%")
RE_DISCARDED    = re.compile(r"Genome discarded:\s*([0-9.]+)%")
RE_TOTAL_MB     = re.compile(r"Total genome length represented:\s*([0-9.]+)\s*Mb")
RE_LINE_BUCKET  = re.compile(r"^(passed|didn't pass)\s*:\s*\d+\s+queries,\s*([0-9.]+)\s*Mb", re.IGNORECASE)

def parse_txt_metrics(path):
    """Return dict with retained_pct, discarded_pct, total_mb, retained_mb, discarded_mb."""
    try:
        with open(path, 'r') as f:
            text = f.read()
    except Exception:
        return None

    retained_pct = None; discarded_pct = None; total_mb = None
    retained_mb = None; discarded_mb = None

    m = RE_RETAINED.search(text)
    if m: retained_pct = float(m.group(1))
    m = RE_DISCARDED.search(text)
    if m: discarded_pct = float(m.group(1))
    m = RE_TOTAL_MB.search(text)
    if m: total_mb = float(m.group(1))

    for status, mb in RE_LINE_BUCKET.findall(text):
        mb_val = float(mb)
        if status.lower().startswith('passed'):
            retained_mb = mb_val
        else:
            discarded_mb = mb_val

    # Infer missing values if needed
    if retained_mb is None and (retained_pct is not None) and (total_mb is not None):
        retained_mb = total_mb * (retained_pct / 100.0)
    if discarded_mb is None and (discarded_pct is not None) and (total_mb is not None):
        discarded_mb = total_mb * (discarded_pct / 100.0)
    if total_mb is None and (retained_mb is not None) and (discarded_mb is not None):
        total_mb = retained_mb + discarded_mb
    if retained_pct is None and (retained_mb is not None) and (total_mb and total_mb > 0):
        retained_pct = 100.0 * retained_mb / total_mb
    if discarded_pct is None and (discarded_mb is not None) and (total_mb and total_mb > 0):
        discarded_pct = 100.0 * discarded_mb / total_mb

    if any(v is None or (isinstance(v, float) and math.isnan(v))
           for v in [retained_pct, discarded_pct, total_mb]):
        return None

    return {
        'retained_pct': retained_pct,
        'discarded_pct': discarded_pct,
        'total_mb': total_mb,
        'retained_mb': retained_mb if retained_mb is not None else total_mb * retained_pct / 100.0,
        'discarded_mb': discarded_mb if discarded_mb is not None else total_mb * discarded_pct / 100.0,
    }

def parse_filename_meta(fname):
    """Extract pair_label, top_mode ('5' or 'all'), and min_coverage (int) from filename."""
    base = os.path.basename(fname)

    # Pair
    pair_label = None
    for prefix, label in PAIR_LABELS.items():
        if base.startswith(prefix):
            pair_label = label
            break

    # top5 / topall
    m = RE_TOP.search(base)
    top_mode = m.group(1).lower() if m else None  # '5' or 'all'

    # cov20 / cov30 / cov50
    m = RE_COV.search(base)
    min_cov = int(m.group(1)) if m else None

    return pair_label, top_mode, min_cov

def main():
    txt_files = sorted(glob.glob(os.path.join(RESULTS_DIR, "*.txt")))
    if not txt_files:
        print(f"[ERROR] No .txt files in {RESULTS_DIR}")
        return

    rows = []
    for path in txt_files:
        pair_label, top_mode, min_cov = parse_filename_meta(path)
        if None in (pair_label, top_mode, min_cov):
            print(f"[SKIP] Unrecognized filename pattern: {os.path.basename(path)}")
            continue
        metrics = parse_txt_metrics(path)
        if metrics is None:
            print(f"[WARN] Could not parse metrics: {path}")
            continue

        rows.append({
            'pair': pair_label,
            'top_mode': 'Top 5' if top_mode == '5' else 'All alignments',
            'min_coverage': min_cov,
            **metrics
        })

    if not rows:
        print("No parsed data. Check filenames and summary format.")
        return

    df = pd.DataFrame(rows)
    df = df.sort_values(['pair', 'min_coverage', 'top_mode'])

    # Save CSV summary
    csv_out = os.path.join(OUT_DIR, "uncut_retained_discarded_summary.csv")
    df.to_csv(csv_out, index=False, float_format='%.2f')
    print(f"Saved summary CSV: {csv_out}")

    # -------- Plot: per pair, retained % vs min_coverage; bars = Top 5 vs All --------
    pairs = list(df['pair'].unique())
    cov_levels = sorted(df['min_coverage'].unique())

    nrows, ncols = len(pairs), 1
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8, 4*nrows), squeeze=False)

    for ax, pair_label in zip(axes[:,0], pairs):
        sub = df[df['pair'] == pair_label]

        top5 = sub[sub['top_mode'] == 'Top 5'].set_index('min_coverage')
        tall = sub[sub['top_mode'] == 'All alignments'].set_index('min_coverage')

        top5_vals = [top5.loc[c, 'retained_pct'] if c in top5.index else float('nan') for c in cov_levels]
        tall_vals = [tall.loc[c, 'retained_pct'] if c in tall.index else float('nan') for c in cov_levels]

        x = range(len(cov_levels))
        width = 0.4

        b1 = ax.bar([i - width/2 for i in x], top5_vals, width, label='Top 5', alpha=0.9)
        b2 = ax.bar([i + width/2 for i in x], tall_vals, width, label='All alignments', alpha=0.9)

        ax.set_xticks(list(x))
        ax.set_xticklabels([f"{c}%" for c in cov_levels])
        ax.set_xlabel('Per-alignment minimum coverage')
        ax.set_ylabel('Genome retained (%)')
        ax.set_title(pair_label)

        # annotate (cleaner placement)
        def annotate(bars):
            for rect in bars:
                h = rect.get_height()
                if math.isnan(h) or h <= 0:
                    continue
                # If bar is tall enough, print inside; otherwise above
                if h > ax.get_ylim()[1] * 0.1:
                    ax.text(
                        rect.get_x() + rect.get_width() / 2, h * 0.5,
                        f"{h:.1f}%", ha='center', va='center',
                        color='black', fontsize=8, fontweight='bold'
                    )
                else:
                    ax.text(
                        rect.get_x() + rect.get_width() / 2, h + (h * 0.05),
                        f"{h:.1f}%", ha='center', va='bottom',
                        fontsize=8, color='black'
                    )

        annotate(b1)
        annotate(b2)

    # One small, shared legend above, top-left, not overlapping
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper left', bbox_to_anchor=(0.01, 1.04), fontsize=8, frameon=False, ncol=1)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plot_out = os.path.join(OUT_DIR, "uncut_retained_by_min_coverage_top5_vs_all.png")
    plt.savefig(plot_out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved plot: {plot_out}")

if __name__ == "__main__":
    main()
