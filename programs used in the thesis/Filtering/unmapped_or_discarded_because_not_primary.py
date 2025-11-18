"""
PAF Alignment Coverage Summation (Total vs. Primary)
====================================================

Description:
    This script processes a Minimap2 PAF alignment file to calculate and compare the total 
    base pairs (bp) aligned across all alignments against the total base pairs aligned 
    in only the **primary** alignments (those marked with the `tp:A:P` tag). 
    This provides a quick measure of the redundancy or fragmentation of alignments, 
    indicating how much aligned sequence is captured by non-primary (supplementary or 
    secondary) hits.

Key Features:
    1. Alignment Type Distinction: Iterates through the PAF and uses the optional Minimap2 
       `tp:A:P` tag to identify primary alignments.
    2. Base Pair Summation: Uses `defaultdict` to efficiently sum the aligned length 
       (`qend - qstart`) for every alignment, grouped by query sequence, in two categories: 
       Total (all hits) and Primary (only 'P' hits).
    3. Quantitative Output: Prints the total sum of aligned Megabases (Mb) for both categories 
       and calculates the **ratio** of primary alignment coverage to total alignment coverage.

Parameters (User Inputs):
    paf (str): Path to the input Minimap2 PAF alignment file.

Outputs:
    - Console Output: Total aligned Mb for all alignments, total aligned Mb for primary 
      alignments, and the percentage ratio of primary aligned Mb to total aligned Mb.
"""
from collections import defaultdict

paf = "/Users/michal/Desktop/thesis_2/data/genomes/unmapped_high_chill_25kb_to_royal_royce.paf"
query_total_aligned = defaultdict(int)
query_primary_aligned = defaultdict(int)

with open(paf) as f:
    for ln in f:
        if not ln.strip() or ln[0] == "#": continue
        p = ln.split("\t")
        qname = p[0]
        qstart, qend = int(p[2]), int(p[3])
        aligned = qend - qstart
        is_primary = any(x == "tp:A:P" for x in p[12:])
        query_total_aligned[qname] += aligned
        if is_primary:
            query_primary_aligned[qname] += aligned

tot_all = sum(query_total_aligned.values())
tot_primary = sum(query_primary_aligned.values())
print(f"Sum of aligned bp across all alignments: {tot_all/1e6:.3f} Mb")
print(f"Sum of aligned bp in primary alignments only: {tot_primary/1e6:.3f} Mb")
print(f"Ratio (primary/all): {100*tot_primary/tot_all:.2f}%")
