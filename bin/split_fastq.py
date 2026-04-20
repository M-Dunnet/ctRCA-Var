#!/usr/bin/env python3

import gzip
import os
import argparse


def open_input(file_path):
    """Open input FASTQ (gz or plain)."""
    if file_path.endswith(".gz"):
        return gzip.open(file_path, "rt")
    return open(file_path, "r")


def open_output(file_path, compress):
    """Open output FASTQ (gz or plain)."""
    if compress:
        return gzip.open(file_path, "wt")
    return open(file_path, "w")


def split_fastq(input_file, reads_per_chunk, output_prefix, compress):
    infile = open_input(input_file)

    chunk_idx = 1
    read_count = 0
    out = None

    base = os.path.basename(input_file)
    if base.endswith(".gz"):
        base = base[:-3]
    base = os.path.splitext(base)[0]

    os.makedirs(output, exist_ok=True)
    
    def new_output():
        nonlocal chunk_idx

        suffix = f".part{chunk_idx:04d}.fastq"
        if compress:
            suffix += ".gz"

        filename = os.path.join(output, base + suffix)
        return open_output(filename, compress)

    out = new_output()

    while True:
        # FASTQ = 4 lines per read
        lines = [infile.readline() for _ in range(4)]

        if not lines[0]:
            break  # EOF

        # Write read
        out.writelines(lines)
        read_count += 1

        # If chunk full → new file
        if read_count >= reads_per_chunk:
            out.close()
            chunk_idx += 1
            read_count = 0
            out = new_output()

    if out:
        out.close()

    infile.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split FASTQ into chunks of N reads")

    parser.add_argument("-i", "--input", required=True, help="Input FASTQ(.gz)")
    parser.add_argument("-n", "--reads", type=int, required=True, help="Reads per chunk")
    parser.add_argument("-o", "--output", required=True, help="Output location")
    parser.add_argument("--gzip", action="store_true", help="Compress output (.gz)")

    args = parser.parse_args()

    split_fastq(
        input_file=args.input,
        reads_per_chunk=args.reads,
        output_prefix=args.output_prefix,
        compress=args.gzip
    )