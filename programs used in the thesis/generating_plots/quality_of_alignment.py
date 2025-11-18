#!/usr/bin/env python3
# qc_minimap2_paf.py
# Analyze Minimap2 PAF alignments: per-hit, per-query, threshold summary, and plots.
# Drop-in: set PAF_PATH, QUERY_FASTA, OUT_DIR below and run in PyCharm.
"""
Minimap2 PAF Best-Hit QC and Visualization Tool
===============================================

Description:
    This script provides a comprehensive quality control analysis of Minimap2 alignments (PAF file). 
    It identifies the single "best" hit for every query sequence and then analyzes the distribution 
    of alignment quality metrics (mismatch rate, MAPQ, coverage) across these best hits. The script 
    generates summary tables and publication-style plots to facilitate the setting of optimal 
    filtering thresholds for genome analysis.

Key Features:
    1. Robust Best-Hit Selection: Determines the single best alignment per query using a hierarchical 
       sorting strategy: **Primary Alignment (tp:A:P) $\to$ Alignment Length (aln\_span) $\to$ 
       Lowest Mismatch Rate $\to$ Highest MAPQ.**
    2. Accurate Mismatch Calculation: Calculates the mismatch rate (%) for all hits, prioritizing 
       information from PAF tags (`de:f` and `NM:i`) for highest accuracy.
    3. Multi-Criteria Gate Analysis: Applies pre-defined quality filters (Loose, Default, Strict) 
       based on MAPQ, Mismatch Rate, Query Coverage, and Alignment Length to quantify the genome 
       Megabases (Mb) retained at each level.
    4. Mismatch Threshold Sweep: Generates a bar plot showing the cumulative volume (Mb and %) of 
       the query genome that would be **discarded** as the maximum allowable mismatch rate is lowered.
    5. Detailed Output: Produces three key tables and multiple plots:
        * `per_hit_metrics.tsv`: All calculated metrics for every alignment in the PAF.
        * `per_query_best.tsv`: Metrics for only the single best alignment of each query.
        * `summary_thresholds.tsv`: Quantitative summary of the impact of the `GATES` criteria.
        * Plots: Histograms of mismatch rate and MAPQ, and a bar chart of total aligned Mb by length bin.

Parameters (User Inputs):
    PAF_PATH (str)       : Path to the input Minimap2 PAF alignment file.
    QUERY_FASTA (str)    : Path to the query sequences (used for accurate query length/coverage).
    OUT_DIR (str)        : Directory for saving all outputs (TSV tables and PNG plots).
    GATES (list)         : Defines the strictness levels (e.g., 'strict', 'default') for the threshold summary.

Outputs:
    - Tables (.tsv): `per_hit_metrics.tsv`, `per_query_best.tsv`, `summary_thresholds.tsv`.
    - Plots (.png): Histograms and comparative bar charts of key quality metrics.
"""

import os, gzip, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -------------------- USER PATHS (EDIT THESE) --------------------
PAF_PATH = "/Users/michal/Desktop/pythonProject1_2/cut_genomes/alignments_with_cut_genomes/alignment_proper_cuts/filtered_alignments/filtered_fasta/hapil_to_florida_core_25kb_29_10.paf"
QUERY_FASTA = "/Users/michal/Desktop/pythonProject1_2/cut_genomes/alignments_with_cut_genomes/alignment_proper_cuts/filtered_alignments/filtered_fasta/filtered_fasta_25kb/hapil_25kb_CORE_only.fasta"
OUT_DIR = "/Users/michal/Desktop/pythonProject1_2/cut_genomes/alignments_with_cut_genomes/alignment_proper_cuts/filtered_alignments/filtered_fasta/quality_alignment_different_sorting_primary_length>mismatch>mapq"
# If reading the FASTA is slow/unreliable (e.g., external/network/iCloud), set False to rely on PAF qlen.
USE_FASTA_LENGTHS = True
# Histogram y-axis tick step for mismatch histogram
HIST_YTICK_STEP = 200

# Thresholds for bar plots (percent mismatch cutoffs)
BAR_THRESHOLDS = [5,10,15,20,25,30,40,50,60]
# -----------------------------------------------------------------

# Gates to summarize pass/fail using MISMATCH RATE (lower is better)
# (name, min_MAPQ, max_mismatch_pct, min_qcov, min_aln_bp)
GATES = [
    ("loose",   20, 5.0,  0.80,  5000),
    ("default", 30, 3.0,  0.90, 10000),
    ("strict",  40, 2.0,  0.95, 15000),
]

def is_gz(p): return str(p).endswith(".gz")
def opengz(p): return gzip.open(p, "rt") if is_gz(p) else open(p, "r")

def fasta_lengths(fa_path):
    lens = {}
    name = None; length = 0
    with opengz(fa_path) as f:
        for ln in f:
            if ln.startswith(">"):
                if name is not None: lens[name] = length
                name = ln[1:].strip().split()[0]; length = 0
            else:
                length += len(ln.strip())
    if name is not None: lens[name] = length
    return lens

def try_load_qlens_from_fasta(path):
    try:
        if not path or not os.path.exists(path):
            print(f"[WARN] FASTA not found: {path!r}. Will use PAF qlen.")
            return {}
        with opengz(path) as f:
            _ = f.readline()
        print(f"[INFO] Reading FASTA lengths from: {path}")
        return fasta_lengths(path)
    except Exception as e:
        print(f"[WARN] Could not read FASTA ({e}). Falling back to PAF qlen.")
        print("[HINT] If this FASTA is on external/network/iCloud storage, copy it locally or set USE_FASTA_LENGTHS=False.")
        return {}

def parse_paf_compute_mismatch(paf_path):
    """
    Parse PAF and compute mismatch rate (%) per hit.
    Priority: de:f -> NM/blen -> (blen - nmatch)/blen
    """
    rows = []
    with opengz(paf_path) as f:
        for ln in f:
            if not ln or ln[0] == "#": continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 12: continue
            try:
                qn, qlen, qs, qe = p[0], int(p[1]), int(p[2]), int(p[3])
                strand = p[4]
                tn, tlen, ts, te = p[5], int(p[6]), int(p[7]), int(p[8])
                nmatch, blen, mapq = int(p[9]), int(p[10]), int(p[11])
            except Exception:
                continue

            nm = None; de = None; cs = None; cg = None
            for opt in p[12:]:
                if   opt.startswith("NM:i:"):
                    try: nm = int(opt.split(":")[-1])
                    except: pass
                elif opt.startswith("de:f:"):
                    try: de = float(opt.split(":")[-1])
                    except: pass
                elif opt.startswith("cs:Z:"):
                    cs = opt[5:]
                elif opt.startswith("cg:Z:"):
                    cg = opt[5:]

            aln_span = max(0, qe - qs)

            # mismatch rate (%)
            if blen > 0:
                if de is not None:
                    mismatch_rate_pct = max(0.0, (1.0 - de) * 100.0)
                elif nm is not None:
                    mismatch_rate_pct = max(0.0, (nm / blen) * 100.0)
                else:
                    mismatch_rate_pct = max(0.0, ((blen - nmatch) / blen) * 100.0)
            else:
                mismatch_rate_pct = np.nan

            qcov = (aln_span / qlen) if qlen > 0 else np.nan

            rows.append({
                "query": qn, "qlen_paf": qlen, "qs": qs, "qe": qe,
                "target": tn, "tlen": tlen, "ts": ts, "te": te,
                "strand": strand,
                "nmatch": nmatch, "blen": blen, "mapq": mapq,
                "aln_span": aln_span, "qcov": qcov,
                "mismatch_rate_pct": mismatch_rate_pct,
                "nm": nm, "de": de, "cs": cs, "cg": cg
            })
    return pd.DataFrame(rows)

def cs_long_nonmatch_stats(cs_str, min_block=10):
    if not cs_str: return 0, 0
    i, n = 0, len(cs_str)
    blocks = []; run_bp = 0
    def flush():
        nonlocal run_bp, blocks
        if run_bp >= min_block: blocks.append(run_bp)
        run_bp = 0
    while i < n:
        c = cs_str[i]
        if c == '=':
            i += 1
            while i < n and cs_str[i].isalpha(): i += 1
            flush()
        elif c in ('*', '+', '-', '~', ':'):
            i += 1; bp = 0
            while i < n and (cs_str[i].isdigit() or cs_str[i].isalpha()):
                if cs_str[i].isdigit():
                    j = i
                    while j < n and cs_str[j].isdigit(): j += 1
                    bp += int(cs_str[i:j]); i = j
                else:
                    bp += 1; i += 1
            run_bp += bp
        else:
            i += 1
    flush()
    return (sum(blocks), len(blocks))

def set_hist_yaxis(ax, counts, step=200, pad=0.04):
    """Auto-cap y-axis at tallest bar, with configurable tick step."""
    ymax = int(np.max(counts)) if len(counts) else 0
    ax.set_ylim(0, ymax if ymax == 0 else ymax * (1.0 + pad))
    if step > 0 and ymax > 0:
        ax.set_yticks(np.arange(0, ymax + step, step))

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # 1) Load PAF
    if not os.path.exists(PAF_PATH):
        raise FileNotFoundError(f"PAF not found: {PAF_PATH}")
    print(f"[INFO] Parsing PAF: {PAF_PATH}")
    df = parse_paf_compute_mismatch(PAF_PATH)
    if df.empty:
        print("[WARN] No alignments parsed. Exiting.")
        return

    # 2) Load FASTA lengths (optional) and finalize qlen/qcov
    qlens = try_load_qlens_from_fasta(QUERY_FASTA) if USE_FASTA_LENGTHS else {}
    df["qlen"] = df["query"].map(qlens).fillna(df["qlen_paf"]).astype(int)
    df["qcov"] = df["aln_span"] / df["qlen"].where(df["qlen"] > 0, pd.NA)

    # 3) Save per-hit metrics
    per_hit_csv = os.path.join(OUT_DIR, "per_hit_metrics.tsv")
    df.to_csv(per_hit_csv, sep="\t", index=False)

    # --- Extract alignment type (tp:A:P/S/other) safely ---
    if "tp" not in df.columns:
        tp_tags = []
        with open(PAF_PATH) as f:
            for ln in f:
                if not ln.strip() or ln[0] == "#":  # skip comments
                    continue
                p = ln.rstrip("\n").split("\t")
                tp_tag = "?"
                for opt in p[12:]:
                    if opt.startswith("tp:A:"):
                        tp_tag = opt.split(":")[-1]
                        break
                tp_tags.append(tp_tag)
        while len(tp_tags) < len(df):
            tp_tags.append("?")
        df["tp"] = tp_tags

    tp_priority = {"P": 0, "S": 1, "?": 2}
    df["tp_rank"] = df["tp"].map(tp_priority).fillna(2)

    # --- Sort per query: primary first → longest → lowest mismatch → highest MAPQ ---
    df_sorted = df.sort_values(
        ["query", "tp_rank", "aln_span", "mismatch_rate_pct", "mapq"],
        ascending=[True, True, False, True, False]
    )
    best = df_sorted.groupby("query", as_index=False).first()

    # 5) Optional mismatch structure via cs:Z
    long_bp, long_cnt = [], []
    for cs in best["cs"]:
        s_bp, s_cnt = cs_long_nonmatch_stats(cs, min_block=10)
        long_bp.append(s_bp); long_cnt.append(s_cnt)
    best["long_nonmatch_bp"] = long_bp
    best["long_block_count"] = long_cnt

    per_query_csv = os.path.join(OUT_DIR, "per_query_best.tsv")
    best.to_csv(per_query_csv, sep="\t", index=False)

    # 6) Threshold sweep summary (using mismatch rate)
    genome_mb = best["qlen"].sum() / 1e6 if len(best) else 0.0
    recs = []
    for name, mq, max_mis, qcv_min, abp_min in GATES:
        mask = (
            (best["mapq"] >= mq) &
            (best["mismatch_rate_pct"] <= max_mis) &
            (best["qcov"] >= qcv_min) &
            (best["aln_span"] >= abp_min)
        )
        n_pass = int(mask.sum())
        mb_queries = best.loc[mask, "qlen"].sum() / 1e6
        mb_aligned = best.loc[mask, "aln_span"].sum() / 1e6
        recs.append({
            "gate": name,
            "N_queries_passing": n_pass,
            "Query_Mb_passing": round(mb_queries, 3),
            "Aligned_Mb_passing": round(mb_aligned, 3),
            "Genome_Mb_total": round(genome_mb, 3),
            "Frac_genome_passing_queries": round((mb_queries/genome_mb) if genome_mb else 0.0, 4),
            "Median_mismatch_pct": round(best.loc[mask, "mismatch_rate_pct"].median() if n_pass else float("nan"), 3),
            "Median_MAPQ": round(best.loc[mask, "mapq"].median() if n_pass else float("nan"), 1),
            "Median_qcov": round(best.loc[mask, "qcov"].median() if n_pass else float("nan"), 4),
            "Median_aln_bp": int(best.loc[mask, "aln_span"].median()) if n_pass else 0,
            "Median_long_nonmatch_bp": int(best.loc[mask, "long_nonmatch_bp"].median()) if n_pass else 0
        })
    pd.DataFrame.from_records(recs).to_csv(
        os.path.join(OUT_DIR, "summary_thresholds.tsv"), sep="\t", index=False
    )

    # 7) Plots

    # 7a) Histogram of mismatch rate (%), with nicer y-axis
    plt.figure()
    counts_mis, bins_mis, _ = plt.hist(best["mismatch_rate_pct"].dropna(), bins=50)
    plt.xlabel("Mismatch rate (%) [best hit per query]")
    plt.ylabel("Count (queries)")
    plt.title("Mismatch rate distribution")
    set_hist_yaxis(plt.gca(), counts_mis, step=HIST_YTICK_STEP, pad=0.04)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "hist_mismatch_rate.png"))
    plt.close()

    # 7b) MAPQ histogram — simple, previous version
    plt.figure()
    best["mapq"].dropna().plot.hist(bins=41)
    plt.xlabel("MAPQ")
    plt.ylabel("Count (queries)")
    plt.title("MAPQ distribution (best hit per query)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "hist_mapq.png"))
    plt.close()

    # 7f) Discarded queries vs mismatch threshold (kept)
    import matplotlib.patheffects as pe

    total_genome_mb = best["qlen"].sum() / 1e6 if len(best) else 0.0
    disc_counts, disc_qmb, disc_pct = [], [], []
    for thr in BAR_THRESHOLDS:
        mask_disc = best["mismatch_rate_pct"] > thr
        disc_counts.append(int(mask_disc.sum()))
        mb_disc = best.loc[mask_disc, "qlen"].sum() / 1e6
        disc_qmb.append(mb_disc)
        disc_pct.append((mb_disc / total_genome_mb * 100.0) if total_genome_mb else 0.0)

    fig, ax = plt.subplots(figsize=(11.5, 6))
    bars = ax.bar([f"{t}%" for t in BAR_THRESHOLDS], disc_counts, width=0.6,
                  color="#1f77b4", edgecolor="black", linewidth=1.1)
    ax.set_xlabel("Mismatch threshold (%)", fontsize=12)
    ax.set_ylabel("Number of discarded queries", fontsize=12)
    ax.set_title("Discarded queries and discarded Mb vs mismatch threshold", fontsize=14, pad=18)

    shift_x = 0.12  # fraction of bar width to shift text right
    for bar, mb, pct in zip(bars, disc_qmb, disc_pct):
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2 + shift_x,
            height + (max(disc_counts) * 0.015),
            f"{mb:.2f} Mb ({pct:.2f}%)",
            ha="center", va="bottom", fontsize=10.5, fontweight="bold",
            color="black", path_effects=[pe.withStroke(linewidth=2, foreground="white")],
        )

    ax.set_xlim(-0.6, len(BAR_THRESHOLDS) - 0.3)
    ax.set_ylim(0, max(disc_counts) * 1.18)
    ax.tick_params(axis="both", labelsize=11)
    ax.grid(axis="y", linestyle="--", linewidth=0.7, alpha=0.6)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "bar_discarded_counts_with_mb_pct.png"), dpi=150)
    plt.close()

    # 7h) Total aligned length (Mb) by alignment-size bins — using ALL alignments (kept)
    import matplotlib.ticker as mticker
    bins_bp = np.array([0, 1_000, 2_000, 5_000, 10_000, 25_000, 50_000,
                        100_000, 250_000, 500_000, np.inf])

    spans = df["aln_span"].dropna().values
    weights_mb = spans / 1e6
    counts_mb, edges = np.histogram(spans, bins=bins_bp, weights=weights_mb)

    def fmt_bp(n):
        if n == np.inf: return "∞"
        if n >= 1_000_000: return f"{n / 1_000_000:.1f}M"
        if n >= 1_000:     return f"{n / 1_000:.0f}k"
        return str(int(n))

    labels = [f"{fmt_bp(edges[i])}-{fmt_bp(edges[i + 1])}" for i in range(len(edges) - 1)]
    labels[-1] = f">{fmt_bp(edges[-2])}"

    fig, ax = plt.subplots(figsize=(11.5, 6))
    bars = ax.bar(np.arange(len(counts_mb)), counts_mb, width=0.8,
                  color="#1f77b4", edgecolor="black", linewidth=1.1)

    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_xlabel("Alignment length bin (bp)")
    ax.set_ylabel("Total aligned length (Mb)")
    ax.set_title("Total aligned length by alignment size (ALL alignments)")

    for b, v in zip(bars, counts_mb):
        if v <= 0: continue
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() * 1.01,
                f"{v:.1f}", ha="center", va="bottom", fontsize=9)

    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda y, _: f"{y:.0f}"))
    ax.margins(y=0.10)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "bar_total_aligned_mb_by_len_bin_all_alignments.png"), dpi=150)
    plt.close()

    print("[OK] Wrote:")
    print(f"  {per_hit_csv}")
    print(f"  {per_query_csv}")
    print(f"  summary_thresholds.tsv + plots in: {OUT_DIR}")

if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)
    main()
