#!/usr/bin/env python3
"""
PAF Mismatch Threshold Sweep Analysis
=====================================

Description:
    This script analyzes a Minimap2 PAF alignment file to determine how many query sequences 
    (and what total length in Megabases) would be discarded at various minimum quality cutoffs.
    The analysis focuses exclusively on the alignment mismatch rate.

Key Features:
    1. Simple Mismatch Calculation: Calculates the mismatch rate (as a percentage) for every 
       alignment based on the standard PAF columns: `num_matches` and `alignment_length`.
    2. Minimum Quality Filter (Permissive): A query sequence is considered 'KEPT' for a given 
       threshold if **at least one** of its alignments has a mismatch rate equal to or below 
       that threshold. This is a highly permissive filtering strategy.
    3. Threshold Sweep: Iterates through a defined list of mismatch percentages to generate 
       a comprehensive summary table.
    4. Quantitative Summary: Reports the number of queries and the total length (Mb and %) 
       of the genome kept and discarded for each threshold.
    5. Visualization: Generates a bar chart showing the number of discarded queries as the 
       mismatch threshold increases, with annotations showing the corresponding total discarded 
       Megabases (Mb) and percentage (%) of the query genome.

Parameters (User Inputs):
    PAF_PATH (str)             : Path to the input Minimap2 PAF alignment file.
    OUT_DIR (str)              : Directory where the summary table and plot will be saved.
    MISMATCH_THRESHOLDS (list) : List of percentage values (e.g., 5, 10, 20) to use as cutoffs 
                                 during the sweep analysis.

Outputs:
    1. CSV: 'mismatch_threshold_summary.csv' (Detailed pass/fail data for each threshold).
    2. Plot: 'discardedQueries_vs_mismatchThreshold_annotatedMBandPCT_shiftClean.png' (Bar chart).
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -------------------- Tunables --------------------
PAF_PATH = '/Users/michal/Desktop/thesis_2/data/genomes/aln_hapil_to_red_gauntlet_asm5.paf'
OUT_DIR  = '/Users/michal/Desktop/pythonProject1_2/plots_coverage_and_mismatch'

# Mismatch thresholds to evaluate (percent)
MISMATCH_THRESHOLDS = [5, 10, 15, 20, 25, 30, 40, 50, 60]
# --------------------------------------------------

os.makedirs(OUT_DIR, exist_ok=True)

# ---- Load PAF (first 12 standard columns, ignore extras) ----
PAF_COLS = [
    'query_name','query_length','query_start','query_end','strand',
    'target_name','target_length','target_start','target_end',
    'num_matches','alignment_length','mapq'
]
df = pd.read_csv(PAF_PATH, sep='\t', header=None, usecols=range(12),
                 names=PAF_COLS, engine='python')

# ---- Mismatch rate (per alignment) ----
df['alignment_length'] = df['alignment_length'].astype('float64')
df['num_matches'] = df['num_matches'].astype('float64')
df['mismatch_rate'] = (1 - (df['num_matches'] / df['alignment_length'])) * 100.0

# ---- Per-query lengths ----
all_queries = df[['query_name','query_length']].drop_duplicates().set_index('query_name')
total_queries = len(all_queries)
total_len_mb  = all_queries['query_length'].sum() / 1_000_000.0

# ---- Sweep thresholds ----
rows = []
g = df.groupby('query_name')['mismatch_rate']
for t in MISMATCH_THRESHOLDS:
    # Query passes if ANY alignment has mismatch <= threshold
    pass_any = g.min() <= t
    kept_q = pass_any[pass_any].index
    disc_q = pass_any[~pass_any].index

    kept_n = len(kept_q)
    disc_n = len(disc_q)

    kept_len_mb = (all_queries.loc[kept_q, 'query_length'].sum() / 1_000_000.0) if kept_n else 0.0
    disc_len_mb = (all_queries.loc[disc_q, 'query_length'].sum() / 1_000_000.0) if disc_n else 0.0

    # Add genome percentages
    kept_pct = (kept_len_mb / total_len_mb) * 100 if total_len_mb > 0 else 0.0
    disc_pct = (disc_len_mb / total_len_mb) * 100 if total_len_mb > 0 else 0.0

    rows.append({
        'mismatch_threshold_pct': t,
        'kept_queries': kept_n,
        'discarded_queries': disc_n,
        'kept_len_mb': round(kept_len_mb, 2),
        'discarded_len_mb': round(disc_len_mb, 2),
        'kept_genome_pct': round(kept_pct, 2),
        'discarded_genome_pct': round(disc_pct, 2)  # ← NEW COLUMN
    })

summary = pd.DataFrame(rows)
summary_path = os.path.join(OUT_DIR, 'mismatch_threshold_summary.csv')
summary.to_csv(summary_path, index=False)
print(f"Saved summary CSV: {summary_path}")


# --------- Plot: discarded queries vs mismatch threshold (annotate discarded Mb + %) ----------
x_labels = [f"{t}%" for t in summary['mismatch_threshold_pct']]
y_counts = summary['discarded_queries'].values
mb_labels = summary['discarded_len_mb'].values
pct_labels = summary['discarded_genome_pct'].values

plt.figure(figsize=(10, 6))
bars = plt.bar(x_labels, y_counts, color="#1f77b4", edgecolor="black", alpha=0.9)
ax = plt.gca()

plt.xlabel('Mismatch threshold (%)')
plt.ylabel('Number of discarded queries')
plt.title('Discarded queries vs mismatch threshold\n(annotation: discarded Mb and % of genome)')
plt.ylim(0, max(y_counts) * 1.15 if len(y_counts) else 1000)

# Annotate with MB discarded and percentage — start at left edge of each bar
if len(y_counts) > 0:
    y_off = max(y_counts) * 0.01
    for rect, mb, pct in zip(bars, mb_labels, pct_labels):
        h = rect.get_height()
        if h > 0:
            x_left = rect.get_x()  # lewa krawędź słupka
            plt.text(
                x_left,               # start napisu przy lewej krawędzi
                h + y_off,            # wysokość jak wcześniej
                f"{mb:.2f} Mb ({pct:.2f}%)",
                ha='left', va='bottom', fontsize=9
            )


plt.tight_layout()
plot_path = os.path.join(OUT_DIR, 'discardedQueries_vs_mismatchThreshold_annotatedMBandPCT_shiftClean.png')
plt.savefig(plot_path, dpi=200)
plt.close()
print(f"Saved plot: {plot_path}")




# ---- Console summary ----
print("\n=== Mismatch threshold sweep (per-alignment mismatch rate) ===")
print(summary.to_string(index=False))
print(f"\nTotal queries considered: {total_queries}")
print(f"Total length represented: {total_len_mb:.2f} Mb")
