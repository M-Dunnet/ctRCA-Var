import pysam
import matplotlib.pyplot as plt
import re
import seaborn as sns
import os
import argparse
import sys


def parse_args():
    """
    Parses command line arguments.
    """
    parser = argparse.ArgumentParser(description='Removes High error reads. Reads are considered if they have more than a set number of misalignments' 
                                     'to the reference genome (default = 3). Deletions or insertions are only conunted as 1, even if multiple bases are affected.',
                                     add_help=True,
                                     prefix_chars='-')
    # Input and Output file locations
    parser.add_argument('input', type=str, action='store',
                        help='The directory location containing input BAM files')
    parser.add_argument('--output', '-o', type=str, action='store', default='',
                        help='The directory location where files will end up. Defaults to the input directory')
    parser.add_argument('--num', '-n', type=int, action='store', default=3,
                        help='Maximum number of misalignments to the reference genome (multi-base deletions or insertions count as 1). Defaults to 3')

    
    return parser.parse_args()


def plot_total_changes(bam_file): ## Currently not used
    """Extract relevant information from a BAM file."""
    samfile = pysam.AlignmentFile(bam_file, "rb")
    changes_total = []

    for read in samfile.fetch():
        read_data = {}
        read_data[read.query_name] = {
            "ref": read.reference_name,
            "start": read.reference_start,
            "end": read.reference_end,
            "cigar": read.cigarstring,
            "MD": read.get_tag('MD'),
            "cs": read.get_tag('cs'),
            "NM": read.get_tag('NM')
        }
        for read, data in read_data.items():
            subs = re.findall(r'\*[acgtn][acgtn]', data['cs'])  # Substitutions; more than 3 subs = go
            ins = re.findall(r'\+[acgtn]+', data['cs'])  # Insertions (any length); more than 1 del or ins = go
            dels = re.findall(r'\-[acgtn]+', data['cs'])  # Deletions (any length); more than 1 del or ins = go;
            total_changes = len(subs) + len(ins) + len(dels)
            if total_changes > 3:
                print(data['cs'])
            changes_total.append(total_changes)
    samfile.close()
    
    data_to_plot = [datapoint if datapoint < 4 else 4 for datapoint in changes_total]
    
    plt.figure(figsize=(8, 5))
    sns.histplot(data_to_plot, discrete=True, stat='percent', )
    plt.title(f"Total Reads = {len(changes_total)}")
    
    plt.show()


def main(args):
    """
    Filters BAM files based on the number of misalignments to the reference genome. Uses the CS tag"
    """
    file_list = [file for file in os.listdir(args.input) if file.endswith('.bam')]

    if len(file_list) == 0:
        print('There are no BAM files in input folder...')
        sys.exit(1)
    print(len(file_list), 'file(s) provided.')

    for file in file_list:
        print(f"Processing {file}")
        output_bam = f"QC_{file}"
        input_file_path = os.path.join(args.input, file)
        output_file_path = os.path.join(args.output, output_bam)
        
        with pysam.AlignmentFile(input_file_path, "rb") as in_bam, \
            pysam.AlignmentFile(output_file_path, "wb", header=in_bam.header) as out_bam:

            for read in in_bam:
                cs_tag = read.get_tag('cs')
                subs = re.findall(r'\*[acgtn][acgtn]', cs_tag)  # Substitutions; more than 3 subs = go
                ins = re.findall(r'\+[acgtn]+', cs_tag)  # Insertions (any length); more than 1 del or ins = go
                dels = re.findall(r'\-[acgtn]+', cs_tag)  # Deletions (any length); more than 1 del or ins = go;
                total_changes = len(subs) + len(ins) + len(dels)
                if total_changes <= 3:
                    out_bam.write(read)

if __name__ == '__main__':
    args = parse_args()
    
    # check input folder
    if not args.input:
        print('Input folder (--input/-i) is required', file=sys.stderr)
        sys.exit(1)
    if not args.input.endswith('/'):
        args.input = args.input + '/'
    

    # check output folder
    if not args.output:
        print('An output location was not provided. Output files will be placed in the input folder')
        args.output = args.input
    if not os.path.exists(args.output):
        os.mkdir(args.output)
    if not args.output.endswith('/'):
        args.output = args.output + '/'

    main(args)
