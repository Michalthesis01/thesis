#!/usr/bin/env python3
"""
FASTA Multi-Scale Splitter
==========================

Description:
    This script processes a directory of FASTA genome files and splits each sequence 
    into smaller, fixed-size chunks. It generates multiple output files for each 
    input genome corresponding to different specified chunk sizes (e.g., 25kb, 50kb).

Key Features:
    1. Batch Processing: Iterates through all .fasta/.fa files in the input directory.
    2. Multi-Scale Output: Generates separate output files for every size in 'chunk_sizes_bp'.
    3. Smart Tail Merging: Prevents the creation of tiny fragments at the end of sequences. 
       If the remainder of a sequence is smaller than a defined fraction of the chunk size, 
       it is merged into the preceding chunk rather than standing alone.

Parameters (User Settings):
    input_dir (str)      : Path to the folder containing source FASTA files.
    out_dir (str)        : Path where processed files will be saved.
    chunk_sizes_bp (list): A list of integers representing chunk lengths in base pairs 
    (e.g., [25000, 50000] for 25kb and 50kb).
    min_tail_frac (float): Threshold (0.0 - 1.0) for merging the final chunk. 
                           Example: 0.30 means if the last chunk is <30% of the target size,
                           it is appended to the previous chunk.
    wrap_width (int)     : Number of characters per line for the output FASTA sequences.

Output:
    Files are named: {original_filename}_{size}kb.fasta
    Headers are renamed: >{original_header}_{chunk_index}_{start_bp}-{end_bp}
"""
import os

# -------- User settings --------
input_dir = '/Users/michal/Desktop/thesis/data/genomes'     # folder with multiple FASTA files
out_dir = '/Users/michal/Desktop/pythonProject1/cut_genomes'
chunk_sizes_bp = [1_000, 25_000, 50_000, 100_000]   # 25 kb, 50 kb, 100 kb
min_tail_frac = 0.30                         # merge tiny tail chunks (<30%)
wrap_width = 60
# --------------------------------


def write_wrapped(fh, seq, width):
    for i in range(0, len(seq), width):
        fh.write(seq[i:i+width] + '\n')


def split_record_to_chunks(header, sequence, chunk_size, min_tail_frac, fh):
    """Split a single FASTA record into chunks, merge short tails."""
    n = len(sequence)
    if n == 0:
        return

    full_chunks = n // chunk_size
    tail = n % chunk_size

    if tail == 0:
        total = full_chunks
    else:
        if tail < int(min_tail_frac * chunk_size) and full_chunks > 0:
            total = full_chunks
        else:
            total = full_chunks + 1

    start = 0
    for idx in range(total):
        if idx < full_chunks - 1:
            end = start + chunk_size
        elif idx == full_chunks - 1:
            if tail != 0 and tail < int(min_tail_frac * chunk_size):
                end = n
            else:
                end = start + chunk_size
        else:
            end = n

        chunk_seq = sequence[start:end]
        out_header = f"{header}_{idx+1}_{start}-{end}"
        fh.write(f">{out_header}\n")
        write_wrapped(fh, chunk_seq, wrap_width)
        start = end


def split_fasta_to_sizes(in_path, out_dir, chunk_sizes, min_tail_frac):
    os.makedirs(out_dir, exist_ok=True)
    base_name = os.path.splitext(os.path.basename(in_path))[0]

    # Prepare writers for each chunk size
    writers = {}
    for cs in chunk_sizes:
        out_path = os.path.join(out_dir, f"{base_name}_{cs//1000}kb.fasta")
        writers[cs] = open(out_path, 'w')

    def flush_record(header, seq):
        if header is None:
            return
        for cs, fh in writers.items():
            split_record_to_chunks(header, seq, cs, min_tail_frac, fh)

    with open(in_path, 'r') as fin:
        current_header = None
        current_seq = []

        for line in fin:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header is not None:
                    flush_record(current_header, ''.join(current_seq))
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        # Flush the last record
        if current_header is not None:
            flush_record(current_header, ''.join(current_seq))

    # Close all writers
    for fh in writers.values():
        fh.close()

    print(f"✅ Finished: {base_name}")
    return [os.path.join(out_dir, f"{base_name}_{cs//1000}kb.fasta") for cs in chunk_sizes]


if __name__ == "__main__":
    fasta_files = [f for f in os.listdir(input_dir)
                   if f.endswith(('.fasta', '.fa'))]

    print(f"Found {len(fasta_files)} FASTA files in {input_dir}")
    for fasta in fasta_files:
        full_path = os.path.join(input_dir, fasta)
        split_fasta_to_sizes(full_path, out_dir, chunk_sizes_bp, min_tail_frac)

    print("\nAll splitting completed! Output FASTAs are in:", out_dir)
