#!/usr/bin/env python3
"""
Alignment Discard Breakdown Visualization Pipeline
==================================================

Description:
    This script processes pre-generated CSV summaries from a complex filtering pipeline (which
    applied minimum coverage, mismatch rate, and total query coverage criteria). It calculates
    the percentage of the query genome discarded at each stage of the filtering process and 
    generates two multi-panel stacked bar charts to visualize these discard reasons.

Key Features:
    1. Multi-Step Discard Analysis: Splits the discarded portion of the genome into three distinct categories:
        * **Min. Coverage Gate:** Queries discarded immediately due to having *zero* alignments passing the minimum coverage threshold.
        * **Mismatch Gate:** Queries discarded because *no* initial alignment passed the minimum mismatch rate threshold.
        * **Total Coverage Fail:** Queries that passed the initial gates but were ultimately discarded because their total merged coverage was insufficient.
    2. Comparative Grid Generation: Produces two 2x4 grid figures, one for the **"Top 5"** alignment strategy and one for **"All alignments"**.
    3. Detailed Visualization: Each panel in the grid represents a specific genome pair and chunk size (e.g., Florida→Royal Royce — 25 kb). The stacked bar chart shows the breakdown of the discarded genome percentage (Mb) across the three reasons, grouped by the minimum coverage level (20%, 30%, 50%).
    4. Data Aggregation: The retained and discarded Megabases are calculated and normalized by the total query length represented in the dataset.

Parameters (User Inputs):
    BASE (str)          : Root directory for the analysis.
    FILTERED_ROOT (str) : Path to the parent folder containing the six subdirectories with filtering results.
    RUN_DIRS (dict)     : Dictionary mapping the 'topN' mode ('top5', 'all') and coverage level (0.20, 0.30, 0.50) to the specific input folder path.
    PATTERNS (list)     : Defines the title and file glob pattern for each panel in the output 2x4 grid.

Outputs:
    - Two PNG Grid Plots: 'grid_discard_breakdown_from_csv_2x4_top5.png' and 
      'grid_discard_breakdown_from_csv_2x4_all.png', visualizing the breakdown of discarded Mb.
"""
import os, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

# -------------------- Paths --------------------
BASE = "/Users/michal/Desktop/pythonProject1_2/cut_genomes/alignments_with_cut_genomes/alignment_proper_cuts"
FILTERED_ROOT = os.path.join(BASE, "filtered_alignments")   # root with the 6 configs
OUT_DIR = os.path.join(BASE, "plots_how_much_is_discarded")
os.makedirs(OUT_DIR, exist_ok=True)

RUN_DIRS = {
    "top5": {
        0.50: f"{FILTERED_ROOT}/filtered_alignments_min.coverage=50_mismatch=30%_topN=5",
        0.30: f"{FILTERED_ROOT}/filtered_alignments_min.coverage=30_mismatch=30%_topN=5",
        0.20: f"{FILTERED_ROOT}/filtered_alignments_min.coverage=20_mismatch=30%_topN=5",
    },
    "all": {
        0.50: f"{FILTERED_ROOT}/filtered_alignments_min.coverage=50_mismatch=30%_topN=all",
        0.30: f"{FILTERED_ROOT}/filtered_alignments_min.coverage=30_mismatch=30%_topN=all",
        0.20: f"{FILTERED_ROOT}/filtered_alignments_min.coverage=20_mismatch=30%_topN=all",
    }
}
# --------------------------------------------------------

# Grid order
PATTERNS = [
    ("Florida→Royal Royce — 1 kb", "*florida*royal*1kb*__coverage_summary_mismatch*.csv"),
    ("Florida→Royal Royce — 25 kb", "*florida*royal*25kb*__coverage_summary_mismatch*.csv"),
    ("Florida→Royal Royce — 50 kb", "*florida*royal*50kb*__coverage_summary_mismatch*.csv"),
    ("Florida→Royal Royce — 100 kb", "*florida*royal*100kb*__coverage_summary_mismatch*.csv"),
    ("Hapil→Red Gauntlet — 1 kb",  "*hapil*gauntlet*1kb*__coverage_summary_mismatch*.csv"),
    ("Hapil→Red Gauntlet — 25 kb",  "*hapil*gauntlet*25kb*__coverage_summary_mismatch*.csv"),
    ("Hapil→Red Gauntlet — 50 kb",  "*hapil*gauntlet*50kb*__coverage_summary_mismatch*.csv"),
    ("Hapil→Red Gauntlet — 100 kb", "*hapil*gauntlet*100kb*__coverage_summary_mismatch*.csv"),
]

COVERAGE_LEVELS = [0.20, 0.30, 0.50]

def find_csv(dir_path, pattern):
    hits = sorted(glob.glob(os.path.join(dir_path, pattern)))
    return hits[0] if hits else None

def load_summary(csv_path):
    usecols = [
        "query_name","query_length","total_query_coverage",
        "coverage_threshold","coverage_pass",
        "mismatch_rate_threshold","alignments_used_for_merging",
        "top_n_considered","passed_both"
    ]
    return pd.read_csv(csv_path, usecols=usecols)

def split_discard(df):
    total_bp = df["query_length"].sum()
    if total_bp <= 0:
        return dict(total_mb=0, keep_pct=0, cov_gate_pct=0, mm_gate_pct=0, cov_fail_pct=0)

    cov_gate_bp = df.loc[df["top_n_considered"] == 0, "query_length"].sum()
    mm_gate_bp  = df.loc[(df["top_n_considered"] > 0) &
                         (df["alignments_used_for_merging"] == 0), "query_length"].sum()
    cov_fail_bp = df.loc[(df["alignments_used_for_merging"] > 0) &
                         (df["coverage_pass"] == "didn't pass"), "query_length"].sum()
    keep_bp     = df.loc[df["passed_both"] == "passed", "query_length"].sum()

    f = lambda bp: (bp / total_bp * 100.0)
    return dict(
        keep_pct     = f(keep_bp),
        cov_gate_pct = f(cov_gate_bp),
        mm_gate_pct  = f(mm_gate_bp),
        cov_fail_pct = f(cov_fail_bp),
        total_mb     = total_bp / 1_000_000.0
    )

def label_seg(ax, x, y0, h, txt):
    if h <= 0: return
    ax.text(
        x, y0 + h/2.0, txt, ha="center", va="center",
        color="white", fontsize=10, fontweight="bold",
        path_effects=[pe.withStroke(linewidth=2.5, foreground="black")]
    )

def build_grid(mode_key, run_dirs):
    fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(20, 9), squeeze=False)
    ymax_seen = 0.0

    for ax, (title, pat) in zip(axes.flatten(), PATTERNS):
        parts = []
        total_mb = None

        for cov in COVERAGE_LEVELS:
            run_dir = run_dirs.get(cov)
            if run_dir is None:
                parts.append(None); continue

            csv_path = find_csv(run_dir, pat)
            if not csv_path:
                parts.append(None); continue

            df = load_summary(csv_path)
            br = split_discard(df)
            parts.append((cov, br))
            total_mb = br["total_mb"]

        if not any(p is not None for p in parts):
            ax.axis("off"); ax.set_title(f"{title}\n(no data)"); continue

        cov_levels = [p[0] for p in parts if p is not None]
        cov_gate   = [p[1]["cov_gate_pct"] for p in parts if p is not None]
        mm_gate    = [p[1]["mm_gate_pct"]  for p in parts if p is not None]
        cov_fail   = [p[1]["cov_fail_pct"] for p in parts if p is not None]
        discarded  = np.array(cov_gate) + np.array(mm_gate) + np.array(cov_fail)

        x = np.arange(len(cov_levels))

        # UPDATED LABELS
        b1 = ax.bar(x, cov_gate, width=0.6, label="Minimum coverage",   color="#8c564b")
        b2 = ax.bar(x, mm_gate,  width=0.6, bottom=cov_gate,
                    label="Mismatch threshold", color="#1f77b4")
        b3 = ax.bar(x, cov_fail, width=0.6,
                    bottom=np.array(cov_gate)+np.array(mm_gate),
                    label="Total coverage", color="#d62728")

        # annotate each segment
        for i in range(len(x)):
            y0 = 0.0
            for h in (cov_gate[i], mm_gate[i], cov_fail[i]):
                if h > 0:
                    label_seg(ax, x[i], y0, h, f"{h:.1f}%")
                y0 += h

        ax.set_xticks(x)
        ax.set_xticklabels([f"{int(c*100)}%" for c in cov_levels])

        # UPDATED X LABEL
        ax.set_xlabel("Per-alignment min coverage")

        ax.set_ylabel("Genome discarded (%)")
        ymax = max(discarded.max()*1.25, 5)
        ymax_seen = max(ymax_seen, ymax)
        tot_txt = f"(total ≈ {total_mb:.1f} Mb)" if total_mb is not None else ""
        ax.set_title(f"{title}\n{tot_txt}")

    for ax in axes.flatten():
        if ax.has_data(): ax.set_ylim(0, ymax_seen)

    handles, labels = axes[0,0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 1.02))

    plt.tight_layout(rect=[0,0,1,0.98])
    out = os.path.join(OUT_DIR, f"grid_discard_breakdown_from_csv_2x4_{mode_key}.png")
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")

def main():
    build_grid("top5", RUN_DIRS["top5"])
    build_grid("all",  RUN_DIRS["all"])

if __name__ == "__main__":
    main()
