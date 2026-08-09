'''
Monte Carlo sampling using Beta-binomial model for INDELS. The in-exact nature of INDELS (i.e. not limited to 1 of 4 possibilies like SNPs) means its very difficult to model 
using a dirichlet distribution. Instead, INDELs are modelled as an Alt (any insertion or deletion) vs Ref (no insertion or deletion) signal with a Beta-binomial.

This class handles those cases
'''

import sys
from collections import defaultdict
import json
import pandas as pd
import numpy as np
from scipy.stats import beta


class VariantAnalyzer:
    '''
    Main class for handling variant data, error data, and computing p-values.
    '''
    
    def __init__(self, test_file, ctrl_file, variant_file, nn=500000):
        '''
        Initilisation of class; requires INS or DEL variant file and reference file. 
        These are the count dictonaries created by the RefSet Class.
        '''
        self.test_file = test_file                                    ## Alpha and beta params for test file
        self.ref_params = self.generate_betabinom_params(ctrl_file)   ## Alpha and beta params for control file
        self.nn = nn                                                  ## Number of monte-carlo iterations
        # self.filtered_variants = self.filter_variants()             ## Filters out known-non-variants

    @staticmethod
    def generate_betabinom_params(data):
        '''
        Generates inital alphas and betas based on position INDEL counts.
        Data is the count dictonary created by the RefSet Class.
        '''
        params = {}
        for position, calls in data.items():
            inital_beta = calls['.'] + 1
            inital_alpha = (sum(calls.values()) -  calls['.']) + 1
            params[position] = {'alpha': inital_alpha, 'beta': inital_beta}
        
        return params

    @staticmethod   ##TODO tidy up the logic....
    def reference_p_value(a_ref, b_ref, n_test, x_test, nn=500000):
        p = np.random.beta(a_ref, b_ref, nn)
        x_sim = np.random.binomial(n_test, p)

        return (np.sum(x_sim >= x_test) + 1) / (nn + 1)


    def run_simulation(self):
        for pos, calls in self.test_file.items():
            test_ref, test_alt = calls['.'], sum(calls.values()) - calls['.']
            n_test = test_ref + test_alt
            ctrl_param = self.ref_params[pos]
            a, b = ctrl_param['alpha'], ctrl_param['beta']
            p = self.reference_p_value(a, b, n_test, test_alt, self.nn)  ## a per-position probability score measuring how consistent the test data is with the reference error model
            if p < 0.01:
                print(pos, test_ref, test_alt, p)
            if pos == 'NC_000017.11_7670685':
                print(pos, test_ref, test_alt, p)
    ## Before simulation is run; check what the variant from del/ins file is. Take any variants that match and run anaylysis.

def _load_reference_file(file_path):
    """
    Load a single reference JSON file and return a dict:
    {position: {base: count, ...}, ...}
    """

    def _collapse_counts(counts):
        """Collapse counts so each base appears once with strand values summed."""
        collapsed = {}
        for key, value in counts.items():
            base = key.split(",")[0].strip("('\" ")
            collapsed[base] = collapsed.get(base, 0) + value
        return collapsed

    with open(file_path, 'r') as f:
        data = json.load(f)

    return {
        position: {
            base: collapsed.get(base, 0)
            for base in collapsed
        }
        for position, counts in data.items()
        for collapsed in [_collapse_counts(counts)]
    }


def filter_test_positions(ctrl_dat, test_dat, variant_dat):
    '''
    Removes positions from the test data that are below a depth of 500, less than 5 alternate counts, or a VAF below 0.1%.
    Furthermore, if positions in the test set are not present in the reference the are excludeded. 
    Finally, Only variants that are in the 'Deletion_variants.tsv' file (pathogenic variants) are retained.
    '''
    def dict_to_df(data):
        '''
        Converts a dictonary to to DF for easy fitlering
        '''
        df = pd.DataFrame.from_dict(data, orient='index').fillna(0)
        df['depth'] = df.sum(axis=1)
        df['alt'] = df['depth'] - df['.']
        df['vaf'] = df['alt'] / df['depth']
        return df
    
    ########################
    ## Convert test and control data to df for easy filtering.
    ########################
    test_df = dict_to_df(test_dat)
    ctrl_df = dict_to_df(ctrl_dat)

    ########################
    ## Mask non-control positions, low depth, and low count positions. 
    ######################## 
    test_df = test_df.loc[test_df.index.intersection(ctrl_df.index)]
    filtered_df = test_df[
        (test_df['depth'] > 500) &
        (test_df['alt'] > 5) &
        (test_df['vaf'] > 0.001)
    ]
    ########################
    ## Filter positions via previous mask
    ######################## 
    filtered_dict = {
        pos: {k: int(v) for k, v in row.items() if v > 0}
        for pos, row in (
            filtered_df
            .drop(columns=['depth', 'alt', 'vaf'])
            .to_dict(orient='index')
            .items()
        )
    }

    ########################
    ## Remove positions no in the target variants list 
    ######################## 
    valid_positions = set(
        variant_dat['Chrom'].str.strip() + '_' + variant_dat['Start.Pos'].astype(str)
    )
    variant_only_dict = {
        k: v for k, v in filtered_dict.items()
        if k in valid_positions
    }

    return filtered_dict


def build_weighted_reference(ctrl_files):
    """
    Builds a depth-weighted, collapsed reference model per position.

    Output format:
    {
        position: {
            '.': weighted_reference_count,
            'alt': weighted_non_reference_count
        }
    }
    """
    ## Container for all control datasets
    ctrl_pool = {}
    ## Load each file and store by index
    for i, file in enumerate(ctrl_files):
        ctrl_dat = _load_reference_file(file)
        ctrl_pool[i] = ctrl_dat

    ## Initialise accumulators
    numerator = defaultdict(lambda: {'.': 0.0, 'alt': 0.0})
    denominator = defaultdict(float)

    ## Loop over all control datasets
    for ctrl_dat in ctrl_pool.values():

        ## Loop over each genomic position
        for pos, counts in ctrl_dat.items():
            ref = counts.get('.', 0)                            ## Extract reference counts safely
            alt = sum(v for k, v in counts.items() if k != '.') ## Collapse all non-reference events into ALT
            depth = ref + alt                                   ## Total depth at this position
            denominator[pos] += depth                           ## Accumulate total depth (for weighting normalisation)
            numerator[pos]['.'] += ref * depth                  ## Weighted accumulation of reference counts
            numerator[pos]['alt'] += alt * depth                ## Weighted accumulation of ALT counts

    ## Final weighted average per position
    weighted_avg = {
        pos: {
            '.': vals['.'] / denominator[pos],
            'alt': vals['alt'] / denominator[pos]
        }
        for pos, vals in numerator.items()
    }

    return weighted_avg


def build_reference_average(ctrl_files):
    """
    Builds a non-weighted (equal contribution) collapsed reference model per position.

    Output format:
    {
        position: {
            '.': mean_reference_count,
            'alt': mean_non_reference_count
        }
    }
    """

    ## Container for all control datasets
    ctrl_pool = {}

    ## Load each file and store by index
    for i, file in enumerate(ctrl_files):
        ctrl_dat = _load_reference_file(file)
        ctrl_pool[i] = ctrl_dat

    ## Initialise accumulators
    numerator = defaultdict(lambda: {'.': 0.0, 'alt': 0.0})
    counts = defaultdict(int)  ## number of datasets contributing per position

    ## Loop over all control datasets
    for ctrl_dat in ctrl_pool.values():

        ## Loop over each genomic position
        for pos, counts_dict in ctrl_dat.items():
            ref = counts_dict.get('.', 0)                            ## reference counts
            alt = sum(v for k, v in counts_dict.items() if k != '.') ## collapse ALT

            numerator[pos]['.'] += ref
            numerator[pos]['alt'] += alt
            counts[pos] += 1

    ## Final average per position
    avg_ref = {
        pos: {
            '.': vals['.'] / counts[pos],
            'alt': vals['alt'] / counts[pos]
        }
        for pos, vals in numerator.items()
    }

    return avg_ref



files = ["/home/dunmi18p/Python_Projects/ctRCA-Var/files/INDELS/Ref_0a_UMI1_pt2_reference_Del_counts.json",
        "/home/dunmi18p/Python_Projects/ctRCA-Var/files/INDELS/Ref_0b_UMI1_pt2_reference_Del_counts.json",
        "/home/dunmi18p/Python_Projects/ctRCA-Var/files/INDELS/Ref_0c_UMI1_pt2_reference_Del_counts.json"]

file2 = "/home/dunmi18p/Python_Projects/ctRCA-Var/_Variant_INDEL_RefSet/QC_ctDNA_ref_025b_Final.bam_UMI1_pt2_reference_Del_counts.json"


ctrl_dat = build_weighted_reference(files)
test_dat = _load_reference_file(file2)
var_dat = pd.read_csv('../files/Target_Mutations/Deletion_variants.csv')

test_dat_clean = filter_test_positions(ctrl_dat, test_dat, var_dat)

obj = VariantAnalyzer(test_dat_clean, ctrl_dat, var_dat)
obj.run_simulation()

## TODO Need list of Ins and Dels from COSMIC.
## TODO Filter list of INS and DELS
## TODO Batch correction like in SNP dirichlet.