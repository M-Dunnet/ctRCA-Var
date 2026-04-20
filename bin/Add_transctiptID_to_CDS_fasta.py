#!/usr/bin/env python3

import argparse

# ----------------------
# Argument parsing
# ----------------------
parser = argparse.ArgumentParser(
    description="Append gene IDs to FASTA headers using CCDS mappings from BED files."
)

parser.add_argument("--canonical_bed", required=True, help="BED file of canonical transcripts")
parser.add_argument("--ccds_bed", required=True, help="BED file containing CCDS mappings")
parser.add_argument("--input_fasta", required=True, help="Input FASTA file")
parser.add_argument("--output_fasta", required=True, help="Output FASTA file")

args = parser.parse_args()

# ----------------------
# Load canonical transcripts
# ----------------------
canonical_transcripts = {
    line.strip().split('\t')[4]
    for line in open(args.canonical_bed)
    if not line.startswith('#')
}

# ----------------------
# Build CCDS → transcript mapping
# ----------------------
ccds_mapping = {}

with open(args.ccds_bed, 'r') as bed_file:
    for line in bed_file:
        if line.startswith('#'):
            continue

        fields = line.strip().split('\t')
        transcript_id = fields[0]
        ccds_identifier = fields[5]

        if transcript_id in canonical_transcripts:
            ccds_mapping[ccds_identifier] = transcript_id

# ----------------------
# Update FASTA headers
# ----------------------
with open(args.input_fasta, 'r') as fasta_file, open(args.output_fasta, 'w') as output_fasta:
    for line in fasta_file:
        if line.startswith('>'):
            fasta_heading = line.strip()
            ccds_identifier = fasta_heading.split(' ')[0].split('_')[2]

            transcript_id = ccds_mapping.get(ccds_identifier)
            if transcript_id:
                fasta_heading = f"{fasta_heading} gene_ID={transcript_id}"

            output_fasta.write(fasta_heading + '\n')
        else:
            output_fasta.write(line)
