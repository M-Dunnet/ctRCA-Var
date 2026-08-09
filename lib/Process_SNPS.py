## Load Imports
import pandas as pd
from pyfaidx import Fasta
## ctRCA module imports
from lib import ctRCA_Funcs as ctRCA
from lib import ReferenceSet
from lib import Alphas
from lib import Dirichlet_Monte_Carlo as dmc

def process_snps(params, file):
    '''
    Function to process SNP data from BAM files.
    '''
    ## Construct refset for the test sample
    print(f"\tLooking for SNP variants")
    snp_refset = ReferenceSet.SNP_ReferenceSet(
        file_path=file,
        bed_path=params['paths']['bedfile'],
        genome_path=params['paths']['genome'],
        collapse_umi=True
    ).reference_set
 
    ## Split the test sample RefSet into positions with potential mutations and positions without. 
    target_positions, non_target_positions = ctRCA.split_refset_dict(snp_refset, params['paths']['snp_targets'])

    ####################################################################################################
    ## First, get target positions for candidate variants and then transform into VCF format 
    ####################################################################################################
    ## Remove strand information from the refset; specifcally in the target positions. 
    target_refset = {
        pos: ctRCA.collapse_strand_counts(counts)
        for pos, counts in target_positions.items()
    }

    ## Convert into a PD.DataFrame; which will become the main VCF file going forward. 
    vcf_df = ctRCA.refset_dict_to_vcf(target_refset, params['paths']['genome'])
    vcf_df = vcf_df[vcf_df["Depth"] >= 500] ## Remove variants with less than 500 total depth
    vcf_df = vcf_df[vcf_df["Alt_Count"] >= 5]   ## Remove varients with less than 5 alternate counts.
    vcf_df = vcf_df[vcf_df["Alt_Prop"] >= 0.001].reset_index(drop=True)    ## Remove variants with a VAF less than 0.001 (0.1%)

    ## Load in BED file and use it to add gene information to the vcf_df
    bedfile = ctRCA.load_bed(params['paths']['bedfile'])
    vcf_df['Gene'] = vcf_df.apply(lambda row: ctRCA.annotate_gene(row, bedfile), axis=1)

    ####################################################################################################
    ## Second, get non-target positions and use these for batch correction and alpha definitions
    ####################################################################################################
    ## Remove strand information from the refset; specifcally in the non-target positions. 
    non_target_refset = {
        pos: ctRCA.collapse_strand_counts(counts)
        for pos, counts in non_target_positions.items()
    }
    ## Define alphas based on control refset data and test data (for batch correction)
    dirichlet_alphas = Alphas.define_alphas(
        test_refset=non_target_refset, 
        genome=Fasta(params['paths']['genome']),
        snp_controls=params['paths']['snp_controls'])  ## Reference Set data for control samples is hardcoded in this script.

    ####################################################################################################
    ## Third, Run Dirichlet Multinomal Monte-Carlo Simulations
    ####################################################################################################
    dmc_analyzer = dmc.VariantAnalyzer(vcf_df, dirichlet_alphas)
    empirical_p_values = dmc_analyzer.run_analysis()
    df = pd.DataFrame(
        empirical_p_values,
        columns=[
            "Chromosome",
            "Position",
            "Type",
            "Ref Allele",
            "Alt Allele",
            "Alt Count",
            "Depth",
            "VAF",
            "P_value", 
            "Info"
        ]
    )

    return df