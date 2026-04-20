canonical_transcripts_bed = '/Users/dunmi18p/Documents/Panel_Checking/ExonLocations414/Coding_exons_Canoinical_Transcipts.bed'
bed_file_path = '/Users/dunmi18p/Desktop/intermedate_files/Transcript_CCDS_ID.bed'
fasta_file_path = '/Users/dunmi18p/Desktop/intermedate_files/Transcript_CDS_sequence.fasta'
output_fasta_path = '/Users/dunmi18p/Desktop/Final_CDS_sequences.fasta'

# Create a dictionary to store mappings between CCDS identifiers and bed file entries
ccds_mapping = {}

# Store only canonical transcripts
canonical_transcripts = {line.strip().split('\t')[4] for line in open(canonical_transcripts_bed) if not line.startswith('#')}
print(canonical_transcripts)

# Read the bed file and store mappings
with open(bed_file_path, 'r') as bed_file:
    for line in bed_file:
        if not line.startswith('#'):
            fields = line.strip().split('\t')
            transcript_id = fields[0]  # Extract the transcript ID
            ccds_identifier = fields[5]  # Extract the CCDS identifier
            if transcript_id not in canonical_transcripts:
                continue
            else:
                ccds_mapping[ccds_identifier] = transcript_id  # Map CCDS identifier to transcript ID

# Modify the fasta file with bed file information
with open(fasta_file_path, 'r') as fasta_file, open(output_fasta_path, 'w') as output_fasta:
    for line in fasta_file:
        if line.startswith('>'):
            fasta_heading = line.strip()
            ccds_identifier = fasta_heading.split(' ')[0].split('_')[2]  # Extract the CCDS identifier
            bed_entry = ccds_mapping.get(ccds_identifier, '')  # Retrieve the corresponding bed file entry
            new_fasta_heading = f'{fasta_heading} gene_ID={bed_entry}' if bed_entry else fasta_heading
            output_fasta.write(new_fasta_heading + '\n')
        else:
            output_fasta.write(line)
