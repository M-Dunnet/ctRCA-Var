
## Load Imports
import argparse
import subprocess
import time
import sys
import os
import pandas as pd
import pathlib as Path

## ctRCA module imports
from lib import ctRCA_Funcs as ctRCA
from lib.Process_INDELS import process_dels, process_ins
from lib.Process_SNPS import process_snps

# from lib import Check_homopolymer_repeats    #TODO needs fixing
# from lib import Exon_status    #TODO needs fixing
# from lib import variant_zygostity    #TODO needs fixing
# from lib import HGVS_changes    #TODO needs fixing
# from lib import Cross_reference_databases    #TODO needs fixing

def parse_args():
    """
    Parses command line arguments.
    Command line arguments must be preceded by its option string (i.e. --input/-i). There are no positional arguments. --input and is required.
    """
    parser = argparse.ArgumentParser(description='Finds mutations from consensus-called amplicon R2C2 reads.',
                                     add_help=True,
                                     prefix_chars='-')
    # Input and Output file locations
    parser.add_argument('--input', '-i', type=str, action='store',
                        help='The directory location containing input BAM files')
    parser.add_argument('--output', '-o', type=str, action='store', default=None,
                        help='The directory location where files will end up. Defaults to the input directory')

    # Optional arguments for additional files, will overwrite the config file paths
    parser.add_argument('--bedfile', '-b', type=str, action='store',
                        help='Path string to the location of a bedfile. Only loci in the bed file will be '
                             'analysed')
    parser.add_argument('--genome', '-g', type=str, action='store',
                        help='Path string to the location of the reference genome.')
    parser.add_argument('--control_refset', '-r', type=str, action='store',
                        help='Path string to the location control sample refset files. This is required to run dirichlet multinomial analysis')
    parser.add_argument('--exon_bed', '-e', type=str, action='store',
                        help='Path string to the location of an exon boundary bedfile. Will determine if variants are in coding regions')
    parser.add_argument('--cds_fasta', '-cds', type=str, action='store',
                        help='Path string to the location of canonical CDS sequences')
    parser.add_argument('--cosmic_variants', '-cv', type=str, action='store',
                        help='Path string to the location of the COSMIC variants filter file (See `Generate_Cosmic_filters` for formatting)')

    # Processing parameters
    parser.add_argument('--UMI_cutoff', '-u', type=int, action='store',
                        help='Sets the cut-off for the minimum number of reads in a UMI family. UMI families with less than this number will be ignored'
                             'Defaults to 1')
    parser.add_argument('--UMI_proportion', '-p', type=float, action='store',
                        help='float value between 0.0 and 1.0 that sets the cut-off for proportion of reads in a UMI '
                             'group containing a mutation. Defaults to 0.66 (2/3 reads need to contain the variant')
    parser.add_argument('--strand_bias_threshold', '-sb', type=float, action='store',
                        help='Strand bias is determined with Fishers exact test. This argument determines the significance threshold.'
                             'defaults to 0.05')
    parser.add_argument('--minimum_coverage', '-c', type=int, action='store',
                        help='Sets the minimum read coverage when determining relevant mutations. Defaults to 50')
    parser.add_argument('--minimum_vaf', '-f', type=float, action='store',
                        help='Sets the minimum frequency when determining relevant mutations. Defaults to 0.001 (0.1%%)')

    return parser.parse_args()


def main(params):
    """
    Main wrapper function that runs the ctRCA pipeline.
    1) Constructs a reference set dictonary for the test sample
    2) Defines Dirichlet multinomial alphas 
    3) Calculates varaint p-values based on emperical dirichlet monte-carlo simulations
    4) Filters variants based on potential pathogenicity using COSMIC and ClinVar databases
    5) Exports data in TSV format
    """

    ## Collect all files in the input dircetory
    file_list = [params['paths']['input']+file for file in os.listdir(params['paths']['input']) if file.endswith('bam')]

    ## Iterate through files // TODO probably should run this in parelle processes...
    for i, file in enumerate(file_list, start=1):
        print(f"Working on file {i} of {len(file_list)}")
        
        ## Process variants of each type
        snp_variants = process_snps(params, file)
        del_variants = process_dels(params, file)
        ins_variants = process_ins(params, file)

        all_variants = pd.concat([snp_variants, del_variants, ins_variants], ignore_index = True)

        basename = os.path.basename(file)
        all_variants.to_csv(f'{params['paths']['output']}{os.path.splitext(basename)[0]}.tsv', sep='\t', index=False)


if __name__ == '__main__':
    ## Start script timer
    start_time = time.time()

    ## Clear terminal window:
    subprocess.run(
        "cls" if os.name == "nt" else "clear",
        shell=True
    )

    ## Parse arguments and config
    args = parse_args()
    config = ctRCA.config_reader()
    final_config = ctRCA.merge_config_with_args(config, args)

    ## Validate required paths
    required_paths = [
        'bedfile',
        'genome',
        'snp_controls',
        'snp_targets',
        'del_controls',
        'del_targets',
        'ins_controls',
        'ins_targets',
        'exon_bed',
        'cds_fasta',
        'cosmic_variants'
    ]
    missing_paths = [p for p in required_paths if not final_config['paths'].get(p)]
    if missing_paths:
        print("ERROR: Missing required paths in config or CLI overrides:", file=sys.stderr)
        for p in missing_paths:
            print(f"  - {p}", file=sys.stderr)
        sys.exit(1)

    ## Ensure input is a directory
    input_dir = final_config['paths']['input']
    if not input_dir.endswith('/'):
        final_config['paths']['input'] += '/'

    ## Default output to input if empty
    if not final_config['paths'].get('output'):
        final_config['paths']['output'] = final_config['paths']['input']

    ## Check the output folder ends with a backslash
    if not final_config['paths']['output'].endswith('/'):
        final_config['paths']['output'] += '/'

    ## Make output dir if it doesnt exist
    if not os.path.exists(config['paths']['output']):
        os.mkdir(config['paths']['output'])

    # Print user settings:
    print(f'---------------------------------------------------------------------------------')
    print(f'RCA Var: Variant caller for Consensus Sequence Rolling Circle Amplification Data')
    print(f'---------------------------------------------------------------------------------\n')
    print(f'Parameters\n-------------------------------------')
    params = final_config['parameters']
    print(f"Minimum coverage\t\t\t{params['minimum_coverage']}")
    print(f"Minimum UMI copies if consensus\t\t{params['UMI_cutoff']}")
    print(f"UMI variant fraction if consensus\t{params['UMI_proportion']}")
    print(f"Minimum VAF\t\t\t\t{params['minimum_vaf']}\n")

    paths = final_config['paths']
    print(f"Input directory\t\t\t\t{paths['input']}")
    print(f"Output directory\t\t\t{paths['output']}")
    print(f"BED file\t\t\t\t{paths.get('bedfile')}")
    print(f"Reference genome\t\t\t{paths.get('genome')}") 
    print(f"Exon BED file\t\t\t\t{paths.get('exon_bed')}")
    print(f"Control sample\t\t\t\t{paths.get('control_refset')}")
    print(f"CDS FASTA\t\t\t\t{paths.get('cds_fasta')}")  
    print(f"COSMIC variants\t\t\t\t{paths.get('cosmic_variants')}") 
    print('\n')

    # Run main
    main(final_config)

    # Print Script Runtime
    end_time = time.time()
    print('\n'+'Runtime:', round(end_time - start_time, 2), 'seconds')