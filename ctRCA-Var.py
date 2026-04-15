
## Load Imports
import argparse
import platform
import time
import sys
import os
import json
import pysam
import numpy as np
import pandas as pd
from pyfaidx import Fasta

## ctRCA module imports
from lib import ctRCA_Funcs as ctRCA
from lib import ReferenceSet
from lib import Alphas
from lib import Dirichlet_Monte_Carlo as dmc
# from lib import Dirichlet_monte_carlo_simulation    #TODO needs fixing
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
    parser.add_argument('--output', '-o', type=str, action='store', default='',
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
    parser.add_argument('--target_variants', '-tv', type=str, action='store',
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
                        help='Sets the minimum frequency when determining relevant mutations. Defaults to 0.001 (0.1%)')

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
    ## Init Genome file once
    genome = Fasta(params['paths']['genome'])
    ## Collect all files in the input dircetory
    file_list = [params['paths']['input']+file for file in os.listdir(params['paths']['input']) if file.endswith('bam')]

    ## Iterate through files // TODO probably should run this in parelle processes...
    for i, file in enumerate(file_list, start=1):
        # print(f"Working on file {i} of {len(file_list)}:")
        # print(file)

        # ## Construct refset for the test sample
        # snp_refset = ReferenceSet.SNP_ReferenceSet(
        #     file_path=file,
        #     bed_path=params['paths']['bedfile'],
        #     genome_path=params['paths']['genome'],
        #     collapse_umi=True
        # ).reference_set

        # #####
        # ## TEMP just saving JSON files
        # file_to_save = file.rsplit('/', 1)[1]
        # with open(f"{file_to_save}_UMI1_pt2_reference_SNP_counts.json", "w") as f:
        #     json.dump(snp_refset, f, indent=4)
        ####
        ## TEMP just using premade JSON files
        with open('/home/dunmi18p/Python_Projects/RCA_Modelling_2/Ref1m_UMI1_reference_SNP_counts.json') as f:
            snp_refset = json.load(f)
        ####

        ## Split the test sample RefSet into positions with potential mutations and positions without. 
        target_positions, non_target_positions = ctRCA.split_refset_dict(snp_refset, "files/Target_Mutations/Sub_Mutations.csv")

        ####################################################################################################
        ## First, get target positions for candidate variants and then transform into VCF format 
        ####################################################################################################
        ## Remove strand information from the refset; specifcally in the target positions. 
        target_refset = {
            pos: ctRCA.collapse_strand_counts(counts)
            for pos, counts in target_positions.items()
        }

        ## Convert into a PD.DataFrame; which will become the main VCF file going forward. 
        vcf_df = ctRCA.refset_dict_to_vcf(target_refset, genome)
        vcf_df = vcf_df[vcf_df["Depth"] >= 500] ## Remove variants with less than 500 total depth
        vcf_df = vcf_df[vcf_df["Alt_Count"] >= 5]   ## Remove varients with less than 5 alternate counts.
        vcf_df = vcf_df[vcf_df["Alt_Prop"] >= 0.001].reset_index(drop=True)    ## Remove variants with a VAF less than 0.001 (0.1%)

        ## Load in BED file and use it to add gene information to the vcf_df
        bedfile = ctRCA.load_bed(params['paths']['bedfile'])
        vcf_df['Gene'] = vcf_df.apply(lambda row: ctRCA.annotate_gene(row, bedfile), axis=1)
        print(vcf_df)
        ####################################################################################################
        ## Second, get non-target positions and use these for batch correction and alpha definitions
        ####################################################################################################
        ## Remove strand information from the refset; specifcally in the non-target positions. 
        non_target_refset = {
            pos: ctRCA.collapse_strand_counts(counts)
            for pos, counts in non_target_positions.items()
        }
        ## Define alphas based on control refset data and test data (for batch correction)
        dirichlet_alphas = Alphas.define_alphas(non_target_refset, Fasta(params['paths']['genome']))  ## Reference Set data for control samples is hardcoded in this script.
  
        ####################################################################################################
        ## Third, Run Dirichlet Multinomal Monte-Carlo Simulations
        ####################################################################################################
        dmc_analyzer = dmc.VariantAnalyzer(vcf_df, dirichlet_alphas)
        empirical_p_values = dmc_analyzer.run_analysis()
        print(empirical_p_values)
        ## Run Dirichlet Multinomial




if __name__ == '__main__':
    ## Start script timer
    start_time = time.time()

    ## Clear terminal window:
    os.system('cls' if platform.system() == 'Windows' else 'clear')

    ## Parse arguments and config
    args = parse_args()
    config = ctRCA.config_reader()
    final_config = ctRCA.merge_config_with_args(config, args)

    ## Validate required paths
    required_paths = [
        'bedfile',
        'genome',
        'control_refset',
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

    ## Print user settings:
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
