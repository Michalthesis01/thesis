#!/usr/bin/env python3
"""
Genome Retention Analysis & Visualization Pipeline
==================================================

Description:
    This script aggregates and visualizes genome retention statistics derived from 
    previous alignment/filtering steps. It parses summary text files to extract 
    "Genome retained" and "Total genome length" metrics, calculates weighted averages, 
    and produces both a summary CSV and a multi-panel comparative plot.

Key Features:
    1. Robust Parsing: extracts retained/discarded percentages and total Mb from 
       unstructured text logs using Regex. Capable of reconstructing missing values 
       if partial data (e.g., only MB but not %) is present.
    2. Weighted Aggregation: Aggregates results across multiple files for the same 
       configuration, weighing percentages by the total genome length represented.
       This ensures that if your data is split into many files of different sizes 
       (e.g., some large chromosomes and some small scaffolds), the large chromosomes 
       dominate the final statistic, which is scientifically correct. 
    3. Parameter Comparison: Groups data by 'min.coverage' (20, 30, 50%) and 
       'topN' mode (Top 5 vs All).
    4. Visualization: Generates a 2x4 grid of bar charts (Rows=Pairs, Cols=Chunk Sizes)
       comparing 'Top 5' vs 'All' filtering strategies with a fixed Y-axis (0-80%).

Parameters (User Inputs):
    BASE_DIR (str)      : Root directory containing the config subfolders.
    CONFIG_FOLDERS (list): List of folder names representing different filter settings
                           (encodes min.coverage and topN settings).
    PAIR_LABELS (dict)  : Mapping of filename prefixes to readable Genome Pair names.
    WANTED_SIZES (list) : Specific chunk sizes (in kb) to filter and plot (e.g., 1, 25, 50, 100).

Outputs:
    1. CSV: 'cut_genomes_retained_discarded_by_pair_and_size.csv' (Aggregated statistics).
    2. Plot: 'retained_by_pair_and_size_top5_vs_all_2x4_y0-80.png' (Bar charts).
"""
import os
import re
import glob
import math
import pandas as pd
import matplotlib.pyplot as plt

# -------------------- Inputs --------------------
BASE_DIR = '/Users/michal/Desktop/pythonProject1_2/cut_genomes/alignments_with_cut_genomes/alignment_proper_cuts/filtered_alignments/'
RESULTS_SUBDIR = 'results'
OUTPUT_DIR = os.path.join(BASE_DIR, 'plots')
os.makedirs(OUTPUT_DIR, exist_ok=True)

CONFIG_FOLDERS = [
    "filtered_alignments_min.coverage=50_mismatch=30%_topN=5",
    "filtered_alignments_min.coverage=30_mismatch=30%_topN=5",
    "filtered_alignments_min.coverage=20_mismatch=30%_topN=5",
    "filtered_alignments_min.coverage=50_mismatch=30%_topN=all",
    "filtered_alignments_min.coverage=30_mismatch=30%_topN=all",
    "filtered_alignments_min.coverage=20_mismatch=30%_topN=all",
]

PAIR_LABELS = {
    'aln_hapil_to_red_gauntlet': 'Hapil→Red Gauntlet',
    'aln_florida_to_royal_royce': 'Florida Brilliance→Royal Royce',
}
PAIR_ORDER = ['Florida Brilliance→Royal Royce', 'Hapil→Red Gauntlet']
WANTED_SIZES = [1, 25, 50, 100]
# ------------------------------------------------

RE_RETAINED = re.compile(r"Genome retained:\s*([0-9.]+)%")
RE_DISCARDED = re.compile(r"Genome discarded:\s*([0-9.]+)%")
RE_TOTAL_MB = re.compile(r"Total genome length represented:\s*([0-9.]+)\s*Mb")
RE_LINE_BUCKET = re.compile(r"^(passed|didn't pass)\s*:\s*\d+\s+queries,\s*([0-9.]+)\s*Mb", re.IGNORECASE)
RE_SIZE = re.compile(r'_(\d+)\s*kb', re.IGNORECASE)

def parse_summary_txt(path):
    retained_pct = None; discarded_pct = None
    total_mb = None; retained_mb = None; discarded_mb = None
    try:
        with open(path, 'r') as f:
            text = f.read()
    except Exception:
        return None

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

def parse_config_meta(folder_name):
    m = re.search(r"min\.coverage=(\d+)", folder_name)
    min_cov = int(m.group(1)) if m else None
    m = re.search(r"topN=(\w+)", folder_name)
    top = m.group(1) if m else None
    return min_cov, top

def get_pair_label(filename):
    base = os.path.basename(filename)
    for prefix, label in PAIR_LABELS.items():
        if base.startswith(prefix):
            return label
    return None

def get_size_kb(filename):
    base = os.path.basename(filename)
    m = RE_SIZE.search(base)
    if m:
        try:
            return int(m.group(1))
        except Exception:
            return 'uncut'
    return 'uncut'

def aggregate_weighted(rows):
    total_mb_sum = sum(r['total_mb'] for r in rows)
    if total_mb_sum <= 0:
        return None, None
    retained_mb_sum = sum(r['retained_mb'] for r in rows)
    discarded_mb_sum = sum(r['discarded_mb'] for r in rows)
    return (100.0 * retained_mb_sum / total_mb_sum,
            100.0 * discarded_mb_sum / total_mb_sum)

def main():
    agg_rows = []

    for cfg in CONFIG_FOLDERS:
        res_dir = os.path.join(BASE_DIR, cfg, RESULTS_SUBDIR)
        if not os.path.isdir(res_dir):
            print(f"[WARN] Missing: {res_dir}")
            continue

        txt_files = sorted(glob.glob(os.path.join(res_dir, "*.txt")))
        if not txt_files:
            print(f"[WARN] No .txt in {res_dir}")
            continue

        buckets = {}
        for path in txt_files:
            pair = get_pair_label(path)
            if pair is None:
                print(f"[SKIP] Unrecognized: {os.path.basename(path)}")
                continue
            size_kb = get_size_kb(path)
            parsed = parse_summary_txt(path)
            if parsed is None:
                print(f"[WARN] Could not parse: {path}")
                continue
            buckets.setdefault((pair, size_kb), []).append(parsed)

        min_cov, top_mode = parse_config_meta(cfg)

        for (pair_label, size_kb), lst in buckets.items():
            retained_w, discarded_w = aggregate_weighted(lst)
            agg_rows.append({
                'pair': pair_label,
                'size_kb': size_kb,
                'config_folder': cfg,
                'min_coverage': min_cov,
                'top_mode': top_mode,
                'retained_pct_weighted': retained_w,
                'discarded_pct_weighted': discarded_w,
                'n_files_aggregated': len(lst),
            })

    if not agg_rows:
        print("No data aggregated.")
        return

    df = pd.DataFrame(agg_rows)
    df = df[df['size_kb'].isin(WANTED_SIZES)]
    df['size_kb'] = pd.Categorical(df['size_kb'], categories=WANTED_SIZES, ordered=True)
    present_pairs = list(df['pair'].dropna().unique())
    pairs = sorted(present_pairs, key=lambda x: PAIR_ORDER.index(x) if x in PAIR_ORDER else 999)
    df = df.sort_values(['pair', 'size_kb', 'min_coverage', 'top_mode'])

    csv_path = os.path.join(OUTPUT_DIR, 'cut_genomes_retained_discarded_by_pair_and_size.csv')
    df.to_csv(csv_path, index=False, float_format='%.2f')
    print(f"Saved summary CSV: {csv_path}")

    # --------- Plot grid: fixed 0–80 y-scale ---------
    sizes = [s for s in WANTED_SIZES if s in set(df['size_kb'])]
    nrows, ncols = len(pairs), len(sizes)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4*ncols + 2, 4*nrows + 1), squeeze=False)
    first_handles_labels = None

    for r, pair_label in enumerate(pairs):
        for c, size in enumerate(sizes):
            ax = axes[r, c]
            sub = df[(df['pair'] == pair_label) & (df['size_kb'] == size)]
            if sub.empty:
                ax.axis('off')
                ax.set_title(f"{pair_label} — {size} kb\n(no data)")
                continue

            cov_levels = sorted(sub['min_coverage'].dropna().unique())
            top5 = sub[sub['top_mode'] == '5'].set_index('min_coverage')
            tall = sub[sub['top_mode'] == 'all'].set_index('min_coverage')
            top5_vals = [float(top5.loc[x, 'retained_pct_weighted']) if x in top5.index else float('nan') for x in cov_levels]
            tall_vals = [float(tall.loc[x, 'retained_pct_weighted']) if x in tall.index else float('nan') for x in cov_levels]

            x = list(range(len(cov_levels)))
            width = 0.4
            b1 = ax.bar([i - width / 2 for i in x], top5_vals, width, label='Top 5', alpha=0.9)
            b2 = ax.bar([i + width / 2 for i in x], tall_vals, width, label='All alignments', alpha=0.9)

            ax.set_xticks(x)
            ax.set_xticklabels([f"{int(v)}%" for v in cov_levels])
            ax.set_xlabel('Per-alignment min coverage')
            ax.set_ylabel('Genome retained (%)')
            ax.set_ylim(0, 85)
            ax.set_yticks(range(0, 81, 10))
            ax.set_title(f"{pair_label}  —  {size} kb")

            def annotate(bars):
                for rect in bars:
                    h = rect.get_height()
                    if math.isnan(h) or h <= 0:
                        continue
                    if h > 8:
                        ax.text(rect.get_x() + rect.get_width() / 2, h * 0.5,
                                f"{h:.1f}%", ha='center', va='center',
                                color='black', fontsize=8, fontweight='bold')
                    else:
                        ax.text(rect.get_x() + rect.get_width() / 2, h + 1.5,
                                f"{h:.1f}%", ha='center', va='bottom',
                                fontsize=8, color='black')

            annotate(b1)
            annotate(b2)

            if first_handles_labels is None:
                first_handles_labels = ax.get_legend_handles_labels()

    if first_handles_labels is not None:
        handles, labels = first_handles_labels
        fig.legend(handles, labels, loc='upper left',
                   bbox_to_anchor=(0.01, 1.06), fontsize=8, frameon=False, ncol=1)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plot_path = os.path.join(OUTPUT_DIR, 'retained_by_pair_and_size_top5_vs_all_2x4_y0-80.png')
    plt.savefig(plot_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved plot: {plot_path}")

if __name__ == "__main__":
    main()
