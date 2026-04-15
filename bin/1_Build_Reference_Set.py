"""
Produces Reference Sets for ctRCA control files.
"""
import os
import json
from lib import Reference_sets
from pyfaidx import Fasta
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="Produces Reference Sets for ctRCA control files."
                                    add_help=True)

    ## Make reference sets arguments
    parser.add_argument('bam', required=True, help='Path to the BAM file')
    parser.add_argument('bed', required=True, help='Path to the BED file') ## This is optional because of BED hardcoding below // will change in the future
    parser.add_argument('genome', required=True, help='Path to the reference genome FASTA file') ## This is optional because of Genome hardcoding below // will change in the future
    parser.add_argument('--exclude', required=False, help='Path to the BED file with positions to exclude')
    parser.add_argument('--depth_threshold', type=int, default=10, help='Depth threshold for variant calling')
    parser.add_argument('--collapse_umi', action='store_true', default=False, help='Whether to collapse UMIs or not')
    parser.add_argument('--output_prefix', default='', help='Prefix for output files')
    
    return parser.parse_args()


def make_reference_sets(args):
    """
    Creates reference sets for SNP counts and saves them as JSON files.
    Relies on `Reference_sets` module for reference set creation.
    """
    bam = args.bam
    exclude = args.exclude
    min_depth = args.depth_threshold
    collapse_umi = args.collapse_umi
    counts_file = (args.output_prefix or "") + "reference_SNP_counts.json"

    ## SNP Counts Reference Set
    counts_ref_set = Reference_sets.SNP_ReferenceSet(bam, bed, genome, depth_threshold=min_depth, collapse_umi=collapse_umi, positions_to_remove=exclude)
    with open(counts_file, "w") as f:
        json.dump(counts_ref_set.reference_set, f, indent=4)


def main():
    args = parse_args()

    ## Clear terminal and print arguments
    os.system('cls' if os.name == 'nt' else 'clear')
    print("=== Running with arguments ===")
    for arg, value in vars(args).items():
        print(f"{arg:20}: {value}")
    print("=" * 35)

    make_reference_sets(args)


if __name__ == "__main__":
    main()

