def remove_duplicates(input_file, output_dir):
    """
    Removes duplicate reads by read ID
    
    :param input_file: String, Input file location
    :param output_dir: String, Output directory location
    :return: None, Writes a new FASTA file
    """
    if not output_dir.endswith('/'):
        output_dir = output_dir+'/'
    
    sequences = {}
    current_id = None
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                current_id = line.strip()[1:]
                if current_id not in sequences:
                    sequences[current_id] = []
            else:
                sequences[current_id].append(line.strip())

    with open(output_dir+'RCA_Var_deduplicated.fasta', 'w') as f:
        for seq_id, seq_data in sequences.items():
            f.write('>' + seq_id + '\n')
            f.write(''.join(seq_data) + '\n')
