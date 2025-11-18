#!/usr/bin/env python3
"""
PAF Alignment Density & Distribution Analysis
=============================================

Description:
    This script analyzes the distribution of previously unmapped sequences (e.g., "failed" reads)
    that have been re-aligned to a target reference genome (PAF format). It calculates alignment 
    density across the genome using fixed-size windows and generates high-quality visualizations 
    to identify "hotspots" where these sequences cluster.

Key Features:
    1. "Best Primary" Filtering: Parses raw PAF alignments to retain only the single best 
       primary alignment per query. It prioritizes alignment length, then mapping quality (MAPQ), 
       and explicitly checks for the Minimap2 primary tag (`tp:A:P`) if available.
    2. Dual-Scale Binning: 
       - Global: 250 kb windows for a concatenated whole-genome view.
       - Local: 1 Mb windows for per-chromosome distributions and violin plots.
    3. Chromosome-Specific Plots:
       - Genome-wide density bar chart (Concatenated).
       - Per-chromosome total alignment bar chart.
       - Violin plots showing the distribution of alignment density per 1 Mb window.
       - Positional histograms for the top K chromosomes with the most aligned bases.

Parameters (User Inputs):
    TARGET_FASTA (str) : Path to the reference genome FASTA (used for chromosome lengths).
    PAF_FAILED_TO_CHR (str): Path to the PAF alignment file (input).
    OUT_DIR (str)      : Directory where plots and TSV tables will be saved.
    CHR_ORDER (str)    : Sorting logic for chromosomes in plots ('fasta', 'length', or 'lex').
    
Outputs:
    - Plots (.png): Density views, violin distributions, and positional histograms.
    - Tables (.tsv): Underlying numerical data for all generated plots.
"""
import os, gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, FuncFormatter

# ================== Inputs ==================
TARGET_FASTA = "/Users/michal/Desktop/thesis_2/data/genomes/farr1.fa"
PAF_FAILED_TO_CHR = "/Users/michal/Desktop/thesis_2/data/genromes/unmapped_high_chill_25kb_to_royal_royce.paf".replace("genromes","genomes")
OUT_DIR = "/Users/michal/Desktop/pythonProject1_2/cut_genomes/alignments_with_cut_genomes/alignment_proper_cuts/filtered_alignments/filtered_fasta/quality_alignment_different_sorting_primary_length_mismatch_mapq/filtered_fasta"
CHR_ORDER = "fasta"   # "length", "lex", or "fasta"
# ============================================

# Binning / plots
DENSITY_BIN_SIZE = 250_000     # concatenated density (250 kb) for genome-wide view
PER_CHR_BIN_SIZE = 1_000_000   # per-chromosome density (1 Mb) for violin & positional plots
FIG_DPI  = 200
TOP_N_PEAKS = 5                 # annotate top-N on violin
TOP_K_CHROMS_POSITIONAL = 5     # how many chromosomes to show positional tracks for

# ---------- helpers ----------
def is_gz(p):
    p = str(p)
    return p.endswith(".gz") or p.endswith(".bgz")

def opengz(p):
    return gzip.open(p, "rt") if is_gz(p) else open(p, "r")

def fasta_lengths_and_order(path):
    lengths = {}
    order = []
    cur = None; L = 0
    with opengz(path) as f:
        for raw in f:
            if raw.startswith(">"):
                if cur is not None:
                    lengths[cur] = L
                cur = raw[1:].split()[0].strip()
                order.append(cur)
                L = 0
            else:
                L += len(raw.strip())
        if cur is not None:
            lengths[cur] = L
    return lengths, order

def only_chromosomes(names):
    out = []
    for n in names:
        low = n.lower()
        if low.startswith("chr") or low.startswith("ch"):
            out.append(n)
    return out

def add_interval_to_bins(start, end, bin_size, acc):
    """Accumulate bp into bins on a single axis; half-open [start, end)."""
    if start is None or end is None or end <= start:
        return
    b0 = start // bin_size
    b1 = (end - 1) // bin_size
    if b0 == b1:
        acc[b0] = acc.get(b0, 0) + (end - start)
        return
    first_end = (b0 + 1) * bin_size
    acc[b0] = acc.get(b0, 0) + (first_end - start)
    for b in range(b0 + 1, b1):
        acc[b] = acc.get(b, 0) + bin_size
    last_start = b1 * bin_size
    acc[b1] = acc.get(b1, 0) + (end - last_start)

def parse_tp_primary(extra_fields):
    """True if tp:A:P present OR tag absent (permissive)."""
    for fld in extra_fields:
        if fld.startswith("tp:A:"):
            return (fld == "tp:A:P")
    return True

# ---------- main ----------
def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # FASTA lengths and sequence sets
    seq_len_all, fasta_order_all = fasta_lengths_and_order(TARGET_FASTA)
    if not seq_len_all:
        print("[ERROR] No sequences parsed from FASTA.")
        return

    chroms = set(only_chromosomes(fasta_order_all))
    contigs = set(seq_len_all.keys()) - chroms

    # Order for chromosomes (plots use only these)
    if CHR_ORDER == "length":
        chrom_order = sorted(chroms, key=lambda c: seq_len_all[c], reverse=True)
    elif CHR_ORDER == "lex":
        chrom_order = sorted(chroms)
    else:
        chrom_order = [c for c in fasta_order_all if c in chroms]

    # Concatenated offsets for chromosomes axis
    offsets = {}
    cum = 0
    for c in chrom_order:
        offsets[c] = cum
        cum += seq_len_all[c]
    genome_span_bp = cum

    # Read PAF and keep best primary per query across ALL targets (chroms + contigs)
    best = {}  # qname -> (tname, tstart, tend, alen, mapq, is_primary)
    seen = kept_any = 0
    with opengz(PAF_FAILED_TO_CHR) as f:
        for ln in f:
            if not ln or ln[0] == "#":
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            seen += 1
            qname = p[0]
            tname = p[5]
            if tname not in seq_len_all:
                continue  # skip if not in FASTA

            try:
                tstart = int(p[7]); tend = int(p[8])
                alen   = int(p[10])
                mapq   = int(p[11])
            except:
                continue

            is_primary = parse_tp_primary(p[12:])
            prev = best.get(qname)
            if prev is None:
                best[qname] = (tname, tstart, tend, alen, mapq, is_primary)
                kept_any += 1
            else:
                _, _, _, alen0, mapq0, is_primary0 = prev
                better = False
                if alen > alen0:
                    better = True
                elif alen == alen0 and mapq > mapq0:
                    better = True
                elif alen == alen0 and mapq == mapq0 and (is_primary and not is_primary0):
                    better = True
                if better:
                    best[qname] = (tname, tstart, tend, alen, mapq, is_primary)
    print(f"[INFO] PAF lines seen: {seen:,}; queries with kept best hit: {kept_any:,}")

    # Accumulate:
    global_bins = {}                            # genome-wide concatenated (chromosomes only, 250 kb)
    per_chr_bp = {c: 0 for c in chrom_order}    # totals
    per_chr_bins = {c: {} for c in chrom_order} # 1 Mb per-chrom bins for violin & positional
    per_contig_bp = {}                          # contigs (TSVs only)
    kept_primary = 0

    for qname, (tname, tstart, tend, alen, mapq, is_primary) in best.items():
        if not is_primary:
            continue
        kept_primary += 1
        span = max(0, tend - tstart)
        if tname in chroms:
            gstart = offsets[tname] + tstart
            gend   = offsets[tname] + tend
            add_interval_to_bins(gstart, gend, DENSITY_BIN_SIZE, global_bins)
            per_chr_bp[tname] += span
            add_interval_to_bins(tstart, tend, PER_CHR_BIN_SIZE, per_chr_bins[tname])
        elif tname in contigs:
            per_contig_bp[tname] = per_contig_bp.get(tname, 0) + span

    print(f"[INFO] Primary alignments kept (after tp filter): {kept_primary:,}")

    # ===== Tables =====
    dens_rows = []
    if global_bins:
        bmin, bmax = min(global_bins.keys()), max(global_bins.keys())
        for b in range(bmin, bmax + 1):
            bp = global_bins.get(b, 0)
            dens_rows.append({
                "bin_index": b,
                "bin_start_bp": b * DENSITY_BIN_SIZE,
                "bin_start_mb": (b * DENSITY_BIN_SIZE) / 1_000_000.0,
                "aligned_from_failed_kb": bp / 1_000.0
            })
    dens_df = pd.DataFrame(dens_rows)

    chr_rows = []
    total_chr_bp = 0
    for c in chrom_order:
        bp = per_chr_bp.get(c, 0)
        total_chr_bp += bp
        chr_rows.append({"chrom": c, "aligned_from_failed_bp": bp, "aligned_from_failed_mb": bp / 1e6})
    chr_df_out = pd.DataFrame(chr_rows)
    total_all_bp = total_chr_bp + sum(per_contig_bp.values())
    chr_df_out["fraction_of_all"] = (chr_df_out["aligned_from_failed_bp"] / max(1, total_all_bp)).fillna(0.0)

    contig_rows = []
    for ctg, bp in sorted(per_contig_bp.items(), key=lambda x: -x[1]):
        contig_rows.append({"contig": ctg, "aligned_from_failed_bp": bp, "aligned_from_failed_mb": bp / 1e6})
    contig_df_out = pd.DataFrame(contig_rows)
    contig_df_out["fraction_of_all"] = (contig_df_out["aligned_from_failed_bp"] / max(1, total_all_bp)).fillna(0.0)

    chr_max_rows = []
    for c in chrom_order:
        vals_mb = [bp/1_000_000.0 for bp in per_chr_bins.get(c, {}).values()]
        c_max = max(vals_mb) if vals_mb else 0.0
        chr_max_rows.append({"chrom": c, "max_mb_per_1Mb_bin": c_max})
    chr_max_df = pd.DataFrame(chr_max_rows)

    # ===== Plots (chromosomes only) =====
    # A) Genome-wide density (250 kb)
    if not dens_df.empty:
        xs_mb = dens_df["bin_start_mb"].to_numpy()
        ys_kb = dens_df["aligned_from_failed_kb"].to_numpy()
        fig, ax = plt.subplots(figsize=(16, 5))
        ax.bar(xs_mb, ys_kb, width=(DENSITY_BIN_SIZE / 1_000_000.0) * 0.98,
               align="edge", edgecolor="none")
        ax.yaxis.set_major_locator(MaxNLocator(nbins=12))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x):,}"))
        ax.grid(axis='y', linestyle='--', linewidth=0.6, alpha=0.6)
        ax.set_xlabel(f"Concatenated chromosome position (Mb) — order: {CHR_ORDER} (chromosomes only)", fontsize=11)
        ax.set_ylabel("Aligned Hapil (kb) per 250 kb window", fontsize=11)
        ax.set_title("Where do the previously-unmapped Hapil sequences land?\nDensity on chromosomes — 250 kb windows",
                     fontsize=12, pad=16)
        n_chr = len(chrom_order)
        if n_chr:
            step = max(1, n_chr // 24)
            x_end_mb = genome_span_bp / 1e6
            chosen_idx = set(range(0, n_chr, step)) | {0, n_chr - 1}
            for i, c in enumerate(chrom_order):
                if i not in chosen_idx: continue
                x = offsets[c] / 1e6
                if 0 < x < (x_end_mb - 0.1):
                    ax.axvline(x=x, color='k', linewidth=0.35, alpha=0.18)
            ax.axvline(x_end_mb - 0.1, color='k', linewidth=0.35, alpha=0.2)
        x_end_mb = genome_span_bp / 1e6
        pad = x_end_mb * 0.005
        ax.set_xlim(-pad, x_end_mb + pad)
        ymax = float(ys_kb.max()) if ys_kb.size else 0.0
        ax.set_ylim(0, (ymax or 1.0) * 1.08)
        plt.tight_layout()
        png_density = os.path.join(OUT_DIR, "failed_to_chr_density_250kb_chr_only.png")
        plt.savefig(png_density, dpi=FIG_DPI, bbox_inches="tight")
        plt.close()

    # B) Horizontal sorted totals (chromosomes only)
    chr_df_sorted = chr_df_out.sort_values("aligned_from_failed_mb", ascending=True)
    fig2, ax2 = plt.subplots(figsize=(12, max(4, 0.28 * len(chr_df_sorted))))
    ax2.barh(chr_df_sorted["chrom"], chr_df_sorted["aligned_from_failed_mb"])
    ax2.set_xlabel("Total aligned Hapil (Mb)")
    ax2.set_ylabel("Chromosome")
    ax2.set_title("Total placement of previously-unmapped Hapil sequence by chromosome (chromosomes only)")
    ax2.grid(axis='x', linestyle='--', linewidth=0.6, alpha=0.6)
    plt.tight_layout()
    png_totals = os.path.join(OUT_DIR, "failed_to_chr_totals_per_chrom_chr_only.png")
    plt.savefig(png_totals, dpi=FIG_DPI, bbox_inches="tight")
    plt.close()

    # C) Violin (1 Mb windows) with Top-N labels in Mb
    data = []
    labels = []
    for c in chrom_order:
        bins_dict = per_chr_bins.get(c, {})
        vals = [bp / 1_000_000.0 for bp in bins_dict.values() if bp > 0]
        data.append(vals if vals else [0.0])
        labels.append(c)
    fig3, ax3 = plt.subplots(figsize=(16, max(5, 0.35 * len(labels))))
    parts = ax3.violinplot(data, vert=False, showmedians=True, widths=0.9)
    for pc in parts['bodies']:
        pc.set_alpha(0.6)
    y_positions = np.arange(1, len(labels) + 1)
    ax3.set_yticks(y_positions)
    ax3.set_yticklabels(labels)
    ax3.set_xlabel("Aligned Hapil (Mb) per 1 Mb window of Royal Royce chromosome")
    ax3.set_ylabel("Chromosome")
    ax3.set_title("Per-chromosome distribution of aligned Hapil per 1 Mb window\n"
                  "(Values > 1 Mb indicate overlapping Hapil chunks)")
    maxes = [max(d) if d else 0.0 for d in data]
    ax3.scatter(maxes, y_positions, marker="|", s=160, linewidths=1.2)
    xmax_data = max(maxes) if maxes else 0.0
    ax3.set_xlim(0, xmax_data * 1.08 + 0.2)
    xlim_right = ax3.get_xlim()[1]
    peak_table = pd.DataFrame({"chrom": labels, "max_mb_per_1Mb_bin": maxes})
    top_peaks = peak_table.sort_values("max_mb_per_1Mb_bin", ascending=False).head(TOP_N_PEAKS)
    for _, row in top_peaks.iterrows():
        chrom = row["chrom"]; mx = float(row["max_mb_per_1Mb_bin"])
        y = labels.index(chrom) + 1
        dx = 0.02 * xlim_right
        text_x = min(mx + dx, 0.97 * xlim_right)
        ax3.annotate(f"{chrom}  max≈{mx:.3f} Mb",
                     xy=(mx, y), xytext=(text_x, y),
                     va="center", fontsize=9, clip_on=True,
                     arrowprops=dict(arrowstyle="-", lw=0.8, alpha=0.7))
    ax3.grid(axis='x', linestyle='--', linewidth=0.6, alpha=0.6)
    plt.tight_layout()
    png_violin = os.path.join(OUT_DIR, "failed_to_chr_bin_density_violin_1Mb_chr_only.png")
    plt.savefig(png_violin, dpi=FIG_DPI, bbox_inches="tight")
    plt.close()

    # D) Combined positional tracks for Top-K chromosomes (one figure; xlabel on EACH subplot)
    topk = chr_df_out.sort_values("aligned_from_failed_mb", ascending=False) \
        .head(TOP_K_CHROMS_POSITIONAL)["chrom"].tolist()

    if topk:
        n = len(topk)
        figH = max(8, 3.6 * n)  # a bit taller so each subplot has room for an xlabel
        fig4, axes = plt.subplots(n, 1, figsize=(18, figH), sharex=False)
        if n == 1:
            axes = [axes]

        for ax, chrom in zip(axes, topk):
            bins = per_chr_bins.get(chrom, {})
            chr_len_bp = seq_len_all[chrom]
            nbins = int(np.ceil(chr_len_bp / PER_CHR_BIN_SIZE))
            xs = np.arange(nbins) * (PER_CHR_BIN_SIZE / 1e6)  # Mb
            ys_mb = np.array([(bins.get(i, 0) / 1e6) for i in range(nbins)])

            ax.bar(xs, ys_mb, width=(PER_CHR_BIN_SIZE / 1e6) * 0.98, align="edge", edgecolor="none")
            ax.set_ylabel("Aligned Mb\nper 1 Mb")
            ax.set_xlabel("Position along chromosome (Mb)")  # <- label on every subplot
            ax.set_title(f"{chrom}: positional density of aligned Hapil (1 Mb windows)")
            ax.grid(axis='y', linestyle='--', linewidth=0.6, alpha=0.6)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=10, prune=None))

            # per-chrom TSV
            per_chrom_tsv = os.path.join(OUT_DIR, f"failed_to_chr_positional_1Mb_{chrom}.tsv")
            pd.DataFrame({"chrom": chrom, "bin_start_mb": xs, "aligned_mb": ys_mb}).to_csv(per_chrom_tsv, sep="\t",
                                                                                           index=False)

        # Tight layout with a little bottom margin so the lowest xlabel isn't clipped
        plt.tight_layout(rect=(0.04, 0.04, 0.98, 0.98))
        png_topk = os.path.join(OUT_DIR, f"failed_to_chr_positional_1Mb_top{n}.png")
        plt.savefig(png_topk, dpi=FIG_DPI, bbox_inches="tight")
        plt.close(fig4)

    # ===== Write TSVs =====
    if not dens_df.empty:
        dens_tsv = os.path.join(OUT_DIR, "failed_to_chr_density_250kb_chr_only.tsv")
        dens_df.to_csv(dens_tsv, sep="\t", index=False)
    chr_tsv = os.path.join(OUT_DIR, "failed_to_chr_totals_per_chrom_chr_only.tsv")
    chr_df_out.to_csv(chr_tsv, sep="\t", index=False)
    chr_max_tsv = os.path.join(OUT_DIR, "failed_to_chr_per_chrom_max_mb_per_1Mb_bin.tsv")
    chr_max_df.to_csv(chr_max_tsv, sep="\t", index=False)
    contig_tsv = os.path.join(OUT_DIR, "failed_to_contig_totals.tsv")
    contig_df_out.to_csv(contig_tsv, sep="\t", index=False)

    # Console summary
    print("[OK] Plots saved:")
    if not dens_df.empty: print("  - failed_to_chr_density_250kb_chr_only.png")
    print("  - failed_to_chr_totals_per_chrom_chr_only.png")
    print("  - failed_to_chr_bin_density_violin_1Mb_chr_only.png")
    if topk:
        print(f"  - failed_to_chr_positional_1Mb_top{len(topk)}.png")
    print("[OK] Tables saved:")
    if not dens_df.empty: print("  - failed_to_chr_density_250kb_chr_only.tsv")
    print("  - failed_to_chr_totals_per_chrom_chr_only.tsv")
    print("  - failed_to_chr_per_chrom_max_mb_per_1Mb_bin.tsv")
    print("  - failed_to_contig_totals.tsv")
    print("  - failed_to_chr_positional_1Mb_{chrom}.tsv for each top chromosome")
    print(f"[INFO] Total Mb on chromosomes: {chr_df_out['aligned_from_failed_mb'].sum():.3f}")
    print(f"[INFO] Total Mb on contigs:     {contig_df_out['aligned_from_failed_mb'].sum() if not contig_df_out.empty else 0.0:.3f}")

if __name__ == "__main__":
    main()
