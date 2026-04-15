#!/bin/bash
# Simple script to run multiple files with minimap2, add an MD-tag, filter by quailty scores (>20), and index the bam files.

# Check the correct number of arguments has been supplied.
if [ $# -ne 4 ]; then
    echo "Usage: $0 <input_directory> <output_directory> <path/to/genome.fa> <path/to/index.mmi>"
    exit 1
fi

# Assign command line arguments to variables
InputDir="$1"
FinalOutputDir="$2"
Genome="$3"
Index="$4"

# Create a temporary output directory
OutputDir="$FinalOutputDir/tmp"
mkdir -p "$OutputDir"

# Iterate through all files with the .fasta extension in the input directory
for file in "$InputDir"/*.fasta; do
    if [ -f "$file" ]; then

        # Generate the output BAM filename
        output_bam="$OutputDir/$(basename "${file%.fasta}.bam")"
        
        # Align the sequences using minimap2 and sort the resulting SAM file using samtools
        ./minimap2 -ax sr --secondary=no --cs "$Index" "$file" |
        samtools sort -o "$output_bam"
        
        # Create a new BAM file with MD tags, and filter by mapping quality (>20)
        output_bam_calmd="$OutputDir/$(basename "${file%.fasta}_calmd.bam")"
        samtools calmd -u "$output_bam" "$Genome" |
        samtools view -q 20 -b -o "$output_bam_calmd"
        
        # Index the BAM file
        samtools index "$output_bam_calmd"
        
        # Move the processed files to the final output directory
        mv -u "$output_bam_calmd" "$FinalOutputDir"
        mv -u "$output_bam_calmd.bai" "$FinalOutputDir"
    fi
done

# Clean up the temporary directory
rm -r "$OutputDir"

# Define the common element to be removed
common_element="_calmd"

# Loop through the files in the source directory and move them to the final output directory
for file in "$FinalOutputDir"/*.bam; do
    if [ -f "$file" ]; then

        # Extract the filename without the path
        filename=$(basename "$file")

        # Remove the common element from the filename
        new_name="${filename//${common_element}/}"

        # Use the 'mv' command to rename the file
        mv -v "$file" "${FinalOutputDir}/${new_name}"
    fi
done

# Repeat the same process for index files (.bai)
for file in "$FinalOutputDir"/*.bai; do
    if [ -f "$file" ]; then

        # Extract the filename without the path
        filename=$(basename "$file")

        # Remove the common element from the filename
        new_name="${filename//${common_element}/}"

        # Use the 'mv' command to rename the file
        mv -v "$file" "${FinalOutputDir}/${new_name}"
    fi
done
