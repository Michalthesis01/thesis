#!/usr/bin/env python3
"""
Genomic Density Plotting of Failed Query Chunks
===============================================

Description:
    This script visualizes the distribution of "failed" query intervals (scaffold chunks) 
    across the entire concatenated reference genome. It aggregates the total length of failed 
    chunks into fixed-size bins along the genome, enabling the identification of regions 
    (scaffolds or local areas) that contain a higher density of low-quality or unmappable 
    sequences.

Key Features:
    1. Scaffold Mapping: Reads the full, **uncut** query FASTA to determine the true length 
       and order of all scaffolds. This establishes a fixed coordinate system for the 
       concatenated genome axis.
    2. Interval Extraction: Reads the output FASTA containing only the **failed** query 
       chunks. It parses the chunk header (e.g., `scaffold_1_1000-26000`) to extract the 
       original scaffold, start, and end coordinates.
    3. Global Binning: Converts the local chunk coordinates into **global genome coordinates** using calculated scaffold offsets. It then accumulates the base pairs of failed intervals 
       into fixed-size bins (`BIN_SIZE`).
    4. Density Visualization: Generates a bar chart where the X-axis represents the concatenated 
       genome (in Mb) and the Y-axis shows the **kilobases (kb) of failed sequence** within 
       each bin.
    5. Scaffold Annotation: Draws faint vertical lines on the plot to delineate the boundaries 
       between major scaffolds, providing crucial context for the observed density peaks.

Parameters (User Inputs):
    UNCUT_FASTA (str)    : Path to the original full-length query genome FASTA.
    FAILED_FASTA (str)   : Path to the FASTA containing only the IDs of failed query chunks 
                           (output from the filter script).
    SCAFF_ORDER (str)    : Determines how scaffolds are arranged on the X-axis ("fasta" default, "length", or "lex").
    BIN_SIZE (int)       : The fixed size (in bp) used for aggregating failed sequence length.

Outputs:
    1. PNG Plot: `failed_density_global_...kb.png` (Visualization of failed sequence density).
    2. TSV Table: `failed_density_global_...kb.tsv` (Tabular data of failed kb per bin).
"""
import os, re, gzip, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, FuncFormatter

# ================== Inputs ==================
# Uncut scaffold FASTA (to compute scaffold lengths & global offsets)
UNCUT_FASTA  = "/Users/michal/Desktop/thesis_2/data/genomes/hapil_20181217_consensus.fa"

# Prefer: FAILED FASTA produced by the filter step (one header per failed query)
FAILED_FASTA = "/Users/michal/Desktop/pythonProject1_2/cut_genomes/alignments_with_cut_genomes/alignment_proper_cuts/filtered_alignments/filtered_fasta/quality_alignment_different_sorting_primary_length_mismatch_mapq/filtered_fasta/hapil_25kb_CORE_only__FAILED_mismatch_gt_60pct_or_unmapped.fasta"
# Optional alternative: plain text with one failed ID per line (leave None to ignore)
FAILED_IDS_TXT = None

OUT_DIR = os.path.dirname(FAILED_FASTA)

# Scaffold ordering for the concatenated axis: "length" (desc) or "lex"
SCAFF_ORDER = "fasta"

# Bin size for global density (bp)
BIN_SIZE = 250_000  # 250 kb (tweak if you want finer/coarser)
FIGSIZE  = (16, 5)
FIG_DPI  = 200
# ============================================

RE_CHUNK = re.compile(r'^(?P<scaf>.+)_(?P<idx>\d+)_(?P<start>\d+)-(?P<end>\d+)$')

def is_gz(p): return str(p).endswith(".gz")
def opengz(p): return gzip.open(p, "rt") if is_gz(p) else open(p, "r")

def scaffold_lengths_from_uncut_fasta(path):
    scaf_len = {}
    cur = None; L = 0
    with opengz(path) as f:
        for raw in f:
            if raw.startswith(">"):
                if cur is not None:
                    scaf_len[cur] = L
                cur = raw[1:].split()[0].strip()
                L = 0
            else:
                L += len(raw.strip())
        if cur is not None:
            scaf_len[cur] = L
    return scaf_len

def read_failed_ids_from_fasta(path):
    ids = []
    with opengz(path) as f:
        for ln in f:
            if ln.startswith(">"):
                ids.append(ln[1:].strip().split()[0])
    return ids

def read_failed_ids_from_txt(path):
    ids = []
    with open(path) as f:
        for ln in f:
            s = ln.strip()
            if s:
                ids.append(s)
    return ids

def parse_chunk_header(h):
    m = RE_CHUNK.match(h.strip())
    if not m:
        return None, None, None
    scaf = m.group("scaf")
    a, b = int(m.group("start")), int(m.group("end"))
    if b < a: a, b = b, a
    return scaf, a, b

def add_interval_to_bins_global(gstart, gend, bin_size, acc):
    """Accumulate bp into global bins (half-open [gstart, gend))."""
    if gstart is None or gend is None or gend <= gstart:
        return
    b0 = gstart // bin_size
    b1 = (gend - 1) // bin_size
    if b0 == b1:
        acc[b0] = acc.get(b0, 0) + (gend - gstart)
        return
    first_end = (b0 + 1) * bin_size
    acc[b0] = acc.get(b0, 0) + (first_end - gstart)
    for b in range(b0 + 1, b1):
        acc[b] = acc.get(b, 0) + bin_size
    last_start = b1 * bin_size
    acc[b1] = acc.get(b1, 0) + (gend - last_start)

def kb_fmt(x, _pos=None): return f"{int(x):,}"

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # 1) Get true scaffold lengths, build scaffold order and global offsets
    scaf_len = scaffold_lengths_from_uncut_fasta(UNCUT_FASTA)
    if not scaf_len:
        print("[ERROR] No scaffolds parsed from uncut FASTA.")
        return
    scaff_df = pd.DataFrame([{"scaffold": s, "length_bp": L} for s, L in scaf_len.items()])

    if SCAFF_ORDER == "length":
        order = scaff_df.sort_values("length_bp", ascending=False)["scaffold"].tolist()
    elif SCAFF_ORDER == "lex":
        order = scaff_df.sort_values("scaffold")["scaffold"].tolist()
    else:
        # Keep the exact order from the uncut FASTA file
        with open(UNCUT_FASTA) as f:
            order = [ln[1:].split()[0].strip() for ln in f if ln.startswith(">")]

    offsets = {}
    cum = 0
    for s in order:
        offsets[s] = cum
        cum += scaf_len[s]
    genome_span_bp = cum

    # 2) Load failed IDs (FAILED FASTA preferred)
    if FAILED_IDS_TXT and os.path.exists(FAILED_IDS_TXT):
        ids = read_failed_ids_from_txt(FAILED_IDS_TXT)
        src = FAILED_IDS_TXT
    else:
        ids = read_failed_ids_from_fasta(FAILED_FASTA)
        src = FAILED_FASTA

    if not ids:
        print("[WARN] No failed IDs found.")
        return

    # 3) Parse intervals and convert to global positions
    total_failed_bp = 0
    bin_acc = {}
    bad = 0
    for rid in ids:
        scaf, a, b = parse_chunk_header(rid)
        if scaf is None or scaf not in scaf_len:
            bad += 1
            continue
        gstart = offsets[scaf] + a
        gend   = offsets[scaf] + b
        total_failed_bp += (b - a)
        add_interval_to_bins_global(gstart, gend, BIN_SIZE, bin_acc)

    if bad:
        print(f"[WARN] Skipped {bad} IDs that didn’t parse or had missing scaffolds.")
    if not bin_acc:
        print("[WARN] No intervals accumulated into bins.")
        return

    # 4) Build table for plotting
    min_bin = min(bin_acc.keys()); max_bin = max(bin_acc.keys())
    rows = []
    for b in range(min_bin, max_bin + 1):
        bp = bin_acc.get(b, 0)
        rows.append({
            "bin_index": b,
            "bin_start_bp": b * BIN_SIZE,
            "bin_start_mb": (b * BIN_SIZE) / 1_000_000.0,
            "failed_kb": bp / 1_000.0
        })
    df = pd.DataFrame(rows)

    # 5) Plot global density with scaffold boundary lines
    xs_mb = df["bin_start_mb"].to_numpy()
    ys_kb = df["failed_kb"].to_numpy()

    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.bar(xs_mb, ys_kb, width=(BIN_SIZE / 1_000_000.0) * 0.98,
           align="edge", color="#1f77b4", edgecolor="none")

    # --- X and Y formatting ---
    ax.yaxis.set_major_locator(MaxNLocator(nbins=12))
    ax.yaxis.set_major_formatter(FuncFormatter(kb_fmt))
    ax.grid(axis='y', linestyle='--', linewidth=0.6, alpha=0.6)

    # --- Axis labels and title ---
    ax.set_xlabel(f"Concatenated genome position (Mb) — scaffolds ordered by {SCAFF_ORDER}", fontsize=11)
    ax.set_ylabel("Failed (kb) per bin", fontsize=11)
    ax.set_title(f"Global density of failed intervals — bin={BIN_SIZE // 1000} kb\nSource: {os.path.basename(src)}",
                 fontsize=12, pad=16)

    # --- Scaffold boundary lines (reduced & drawn inside the frame) ---
    n_scaff = len(order)
    if n_scaff:
        # draw every ~step-th boundary; adaptive so we cap at ~100 lines total
        step = max(10, n_scaff // 100)  # at most ~100 lines; never denser than every 10th
        x_end_mb = genome_span_bp / 1e6

        # always include first and last scaffold boundaries
        chosen_idx = set(range(0, n_scaff, step)) | {0, n_scaff - 1}
        for i, s in enumerate(order):
            if i not in chosen_idx:
                continue
            x = offsets[s] / 1e6
            if 0 < x < (x_end_mb - 0.1):  # keep slightly inside the right edge
                ax.axvline(x=x, color='k', linewidth=0.3, alpha=0.12)

        # mark end slightly inside to avoid the dark band on the edge
        ax.axvline(x_end_mb - 0.1, color='k', linewidth=0.3, alpha=0.15)

    # --- Add small (balanced) horizontal margins ---
    x_end_mb = genome_span_bp / 1e6
    pad = x_end_mb * 0.005  # 0.5% of genome span
    ax.set_xlim(-pad, x_end_mb + pad)  # symmetrical small gap left/right

    ymax = float(ys_kb.max()) if ys_kb.size else 0.0
    ax.set_ylim(0, (ymax or 1.0) * 1.08)

    # ... plotting code above ...

    plt.tight_layout()
    png_path = os.path.join(OUT_DIR, f"failed_density_global_{BIN_SIZE // 1000}kb.png")
    plt.savefig(png_path, dpi=FIG_DPI, bbox_inches="tight")
    plt.close()

    tsv_path = os.path.join(OUT_DIR, f"failed_density_global_{BIN_SIZE // 1000}kb.tsv")
    df.to_csv(tsv_path, sep="\t", index=False)

    print(f"[OK] Saved plot  → {png_path}")
    print(f"[OK] Saved table → {tsv_path}")
    print(f"[INFO] Total failed interval sum: {total_failed_bp / 1e6:.3f} Mb")
    print(f"[INFO] Genome span (concatenated): {genome_span_bp / 1e6:.2f} Mb")


if __name__ == "__main__":
    main()
