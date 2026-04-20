#!/bin/bash

# Usage: ./trim_adapters.sh <input_dir> <output_dir>

set -e

# Check arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"

# Check input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# Create output directory if needed
mkdir -p "$OUTPUT_DIR"

echo "Input directory:  $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Loop through all .fasta files
for file in "$INPUT_DIR"/*.fasta; do
    # Skip if no files match
    [ -e "$file" ] || continue

    sample=$(basename "$file" .fasta)

    trimmed="$OUTPUT_DIR/${sample}.trimmed.fasta"
    untrimmed="$OUTPUT_DIR/${sample}.untrimmed.fasta"

    echo "Processing $sample ..."

    cutadapt \
        -a GCATTCGAGTCAT...CAAAACGCAATACTGTACTGGAGATCGGA \
        --untrimmed-output "$untrimmed" \
        -o "$trimmed" \
        "$file"
done

echo ""
echo "All samples processed!"
