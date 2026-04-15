#!/bin/bash

# Directory containing your FASTA files
INPUT_DIR="/projects/health_sciences/bms/biochemistry/guilford_group/MichaelDunnet/RCA_cfDNA_standards/Alignment_Metrics_UpdatedStrandInfo/Repeat-Level-Analysis/"
OUTPUT_DIR="/projects/health_sciences/bms/biochemistry/guilford_group/MichaelDunnet/RCA_cfDNA_standards/Alignment_Metrics_UpdatedStrandInfo/Repeat-Level-Analysis/Adapter_Trim/"
mkdir -p "$OUTPUT_DIR"

# Loop through all .fasta files
for file in "$INPUT_DIR"/*.fasta; do
    # Get sample name without extension
    sample=$(basename "$file" .fasta)

    # Define outputs
    trimmed="$OUTPUT_DIR/${sample}.trimmed.fasta"
    untrimmed="$OUTPUT_DIR/${sample}.untrimmed.fasta"

    echo "Processing $sample ..."

    cutadapt \
        -a GCATTCGAGTCAT...CAAAACGCAATACTGTACTGGAGATCGGA \
        --untrimmed-output "$untrimmed" \
        -o "$trimmed" \
        "$file"
done

echo "All samples processed!"
