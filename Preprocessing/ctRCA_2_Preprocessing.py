import argparse
import os
import sys
import platform
from preprocessing_lib import Define_UMI, Remove_duplicates


def parse_args():
    """
    Parses command line arguments.

    Command line arguments must be preceded by its option string (i.e. --input/-i). There are no positional arguments, but --input and --bedfile are
    required.

    Returns:
        Parsed commandline arguments
    """
    parser = argparse.ArgumentParser(description='De-duplicates and defines read UMI.',
                                     add_help=True,
                                     prefix_chars='-')
    # Input and Output file locations
    parser.add_argument('input', type=str, action='store',
                        help='The directory location containing input Fasta files')
    parser.add_argument('output', type=str, action='store', default='',
                        help='The directory location where files will end up. Defaults to the input directory')

    return parser.parse_args()


def main(arg):
    # Create list to hold input files
    file_list = []
    
    # Make output folder if it doesn't exist
    if not os.path.exists(arg.output):
        os.mkdir(arg.output)
    
    # Add files to file list
    for file in os.listdir(arg.input):
        if file.endswith('.fasta') or file.endswith('fa'):
            file_list.append(file)
    
    # Check file list contains FASTA files
    if len(file_list) == 0:
        print('There are no Fasta files in input folder...')
        sys.exit(1)
    
    for file in file_list:
        print(f"Working on file: {file}")
        filename = os.path.basename(file)
        Remove_duplicates.remove_duplicates(arg.input + file, arg.output)
        Define_UMI.define_umi(arg.output + 'RCA_Var_deduplicated.fasta', arg.output, filename)
        
        os.remove(arg.output + 'RCA_Var_deduplicated.fasta')
        os.remove(arg.output + 'RCA_Var_deduplicated.fasta.fai')
        

if __name__ == '__main__':
    
    system_platform = platform.system()
    if system_platform == 'Windows':
        os.system('cls')
    else:
        os.system('clear')
    
    # Parse arguments
    args = parse_args()
    
    # Check an input directory was provided
    if not args.input:
        print('Input folder is required', file=sys.stderr)
        sys.exit(1)
    
    # Check the input folder ends with a backslash
    if not args.input.endswith('/'):
        args.input = args.input + '/'
    
    # If an output folder is not supplied, make the input folder the output
    if not args.output:
        print('An output location was not provided. Output files will be placed in the input folder')
        args.output = args.input
    
    # Check the output folder ends with a backslash
    if not args.output.endswith('/'):
        args.output = args.output + '/'
    
    main(args)
    