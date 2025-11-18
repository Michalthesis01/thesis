
#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------- Paths ----------
PAF_PATH = '/Users/michal/Desktop/thesis_2/data/genomes/aln_hapil_to_red_gauntlet_asm5.paf'
OUT_DIR  = '/Users/michal/Desktop/pythonProject1_2/plots_coverage_and_mismatch'

# ---------- Tunables ----------
PER_ALIGN_COV_MIN = 10.0      # percent; skip tiny alignments
HIST_BIN_WIDTH    = 5         # histogram bin width in percentage points

os.makedirs(OUT_DIR, exist_ok=True)

# ---------- Helpers ----------
def merge_intervals(intervals):
    """Merge [start, end) intervals respecting PAF end-exclusive coords."""
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: x[0])
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        # end-exclusive merge (no +1)
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged

def total_coverage_pct_for_query(sub):
    """Compute unique covered % of the query after merging [start,end) blocks."""
    qlen = int(sub['query_length'].iloc[0])
    if qlen <= 0:
        return 0.0
    intervals = list(zip(sub['query_start'].astype(int), sub['query_end'].astype(int)))
    merged = merge_intervals(intervals)
    uniq_bp = sum(e - s for s, e in merged) if merged else 0
    return (uniq_bp / qlen) * 100.0

# ---------- Load PAF robustly (first 12 core columns only) ----------
PAF_COLS = [
    'query_name','query_length','query_start','query_end','strand',
    'target_name','target_length','target_start','target_end',
    'num_matches','alignment_length','mapq'
]
df = pd.read_csv(PAF_PATH, sep='\t', header=None, usecols=range(12), names=PAF_COLS, engine='python')

# ---------- Optional per-alignment coverage gate (no mismatch filtering) ----------
span_len = (df['query_end'] - df['query_start']).astype('int64')
span_len = span_len.where(span_len > 0, other=1)     # avoid zero division later
df['query_coverage'] = (span_len / df['query_length']) * 100.0
df_filt = df.loc[df['query_coverage'] > PER_ALIGN_COV_MIN].copy()

# ---------- Per-query total coverage ----------
tot_cov = (
    df_filt.groupby('query_name', as_index=False)
           .apply(lambda g: total_coverage_pct_for_query(g), include_groups=False)
           .rename(columns={None: 'total_coverage_pct'})
)

# Keep one query_length per query for MB annotations
qlens = df_filt.groupby('query_name', as_index=False)['query_length'].first()
tot_cov = tot_cov.merge(qlens, on='query_name', how='left')

# ---------- Histogram + annotate top-3 bins by count (tie-break by MB) ----------
bins = np.arange(0, 100 + HIST_BIN_WIDTH, HIST_BIN_WIDTH)
counts, edges = np.histogram(tot_cov['total_coverage_pct'], bins=bins)
centers = (edges[:-1] + edges[1:]) / 2

# Map each query to a bin index
bin_idx = np.digitize(tot_cov['total_coverage_pct'], bins=edges, right=False) - 1
bin_idx = np.clip(bin_idx, 0, len(edges) - 2)

# Sum query lengths per bin (bp -> MB)
sum_bp_per_bin = np.zeros(len(edges) - 1, dtype=float)
for i, qlen_bp in zip(bin_idx, tot_cov['query_length'].to_numpy()):
    sum_bp_per_bin[i] += float(qlen_bp)
sum_mb_per_bin = sum_bp_per_bin / 1e6

# Choose top 3 bins by count, breaking ties by MB (desc)
# Use lexsort with primary key counts, secondary key MB
order = np.argsort(np.lexsort((-sum_mb_per_bin, counts)))[::-1]
top3 = [idx for idx in order if counts[idx] > 0][:3]
top3 = set(top3)

# ---------- Plot ----------
plt.figure(figsize=(10, 6))
bars = plt.bar(centers, counts, width=np.diff(edges), align='center')
plt.xlabel('Total coverage per query (%)')
plt.ylabel('Number of queries')
plt.title('Distribution of total coverage per query')
plt.xticks(bins)

# Annotate only top 3 bins with total query size in MB
y_offset = max(counts) * 0.01 if counts.max() > 0 else 1
for idx, (rect, c) in enumerate(zip(bars, counts)):
    if idx in top3:
        mb = sum_mb_per_bin[idx]
        plt.text(rect.get_x() + rect.get_width()/2,
                 rect.get_height() + y_offset,
                 f'{mb:.1f} MB',
                 ha='center', va='bottom', fontsize=8)

plt.tight_layout()
out_path = os.path.join(OUT_DIR, 'total_coverage_per_query_top3MB.png')
plt.savefig(out_path, dpi=200)
plt.close()

print(f"Saved Plot: {out_path}")
print(f"Alignments kept: {len(df_filt):,} / {len(df):,}")
print(f"Queries represented: {len(tot_cov):,}")
