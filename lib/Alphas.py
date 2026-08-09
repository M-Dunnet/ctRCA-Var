import json
import numpy as np
from pyfaidx import Fasta
from collections import defaultdict
from scipy.special import digamma, polygamma
import matplotlib.pyplot as plt
import seaborn as sns


#########################
## GLOBALS 
#########################
VALID_BASES = ['A', 'C', 'G', 'T']
#########################

def _collapse_counts(counts):
    """Collapse counts so each base appears once with strand values summed."""
    collapsed = {}
    for key, value in counts.items():
        base = key.split(",")[0].strip("('\" ")
        collapsed[base] = collapsed.get(base, 0) + value
    return collapsed



def _load_reference_file(file_path):
    """
    Load a single reference JSON file and return a dict:
    {position: {base: count, ...}, ...}
    """
    with open(file_path, 'r') as f:
        data = json.load(f)

    return {
        position: {
            base: collapsed.get(base, 0)
            for base in VALID_BASES
        }
        for position, counts in data.items()
        for collapsed in [_collapse_counts(counts)]
    }


def _estimate_alphas(VALID_BASES, FILE_PATHS, genome):
    # 1. Pre-load all data into a 3D NumPy array: [File, Position, Base]
    # This avoids re-reading files for every locus.
    
    # First, get a stable list of all loci present in the first file
    first_ref = _load_reference_file(FILE_PATHS[0])
    loci = list(first_ref.keys())
    locus_to_idx = {locus: i for i, locus in enumerate(loci)}
    
    num_files = len(FILE_PATHS)
    num_loci = len(loci)
    num_bases = len(VALID_BASES)
    
    # Create a large zeroed array (Memory check: ~40MB per 1M loci)
    # Shape: (Files, Loci, Bases)
    counts_matrix = np.zeros((num_files, num_loci, num_bases))

    for f_idx, path in enumerate(FILE_PATHS):
        ref_data = _load_reference_file(path)
        for l_idx, locus in enumerate(loci):
            pos_counts = ref_data.get(locus, {})
            # Add + 1e-12 pseudocount here to match your 'build_initial_alphas' logic
            counts_matrix[f_idx, l_idx, :] = [pos_counts.get(b, 0) + 1e-12 for b in VALID_BASES]    #1e-12
    
    
    # # 2. Vectorized Math (Replaces build_initial_alphas)
    # # Kappa is the sum across bases: Shape (Files, Loci)
    kappas = counts_matrix.sum(axis=2) 
    median_depth = np.median([d for d in kappas.flatten() if d > 1000])

    depths = counts_matrix.sum(axis=2, keepdims=True)
    proportions = counts_matrix / depths  # (Files, Loci, Bases)
    mu = np.mean(proportions, axis=0)  # (Loci, Bases)
    var = np.var(proportions, axis=0, ddof=1)  # (Loci, Bases)
    # avoid division by zero / numerical issues
    eps = 1e-12
    var = np.clip(var, eps, None)
    alpha0_per_base = (mu * (1 - mu)) / var - 1
    alpha0 = np.median(alpha0_per_base, axis=1)  # (Loci,)
    init_alphas = mu * alpha0[:, np.newaxis]  # (Loci, Bases)
    init_alphas = np.clip(init_alphas, 1e-8, None)

    # 3. Aggregate by Reference Base
    global_alpha_p = defaultdict(list)
    base_to_idx = {base: i for i, base in enumerate(VALID_BASES)}

    for l_idx, locus in enumerate(loci):
        chromosome, pos_str = locus.rsplit('_', 1)
        ref_base = str(genome[chromosome][int(pos_str) - 1]).upper()
        
        # Calculate alpha_p (normalized alphas for this position)
        alphas = init_alphas[l_idx]
        total = alphas.sum()
        if total < 500: continue    ## Use only positions with more than 500 depth
        alpha_p = alphas / total

        ## Filter Out High-Error Loci ---
        ## Get the index for the reference base (e.g., 'A' -> 0)
        ref_idx = base_to_idx.get(ref_base)
        
        ## Identify "other" indices (everything except the ref_base)
        other_indices = [i for i in range(len(VALID_BASES)) if i != ref_idx]
        
        ## Check if any proportion in the "other" bases is > 0.01
        ## .any() returns True if at least one mismatch base is > 1%
        if (alpha_p[other_indices] > 0.005).any():   
            continue
        
        ## Convert back to dict for existing downstream logic
        alpha_p_dict = dict(zip(VALID_BASES, alpha_p))
        global_alpha_p[ref_base].append((alpha_p_dict))

    gbl_alpha = defaultdict(list)
    for base, alphas_p in global_alpha_p.items():

        ## Calculate means and variances
        n = len(alphas_p)
        keys = alphas_p[0].keys()
        
        means = {
            k: sum(a[k] for a in alphas_p) / n 
            for k in keys
        }

        variances = {
            k: sum((a[k] - means[k])**2 for a in alphas_p) / (n-1)
            for k in keys
        }
        precision_params = {
            k: (means[k] * (1 - means[k]) / variances[k]) - 1
            if variances[k] > 0 else 0  # Avoid division by zero
            for k in keys
        }
        
        mean_kappa = np.mean(list(precision_params.values())) ## Going to take mean because it is more conservative...
        global_alpha = [mean_kappa * p for p in means.values()]
        gbl_alpha[base] = global_alpha

    return gbl_alpha, init_alphas, median_depth


def _calculate_test_errs(test_dict, genome):
    """
    Gets mean error rate from the test set. 
    """
    items = list(test_dict.items())  ## list of positions, index-aligned with counts
    loc = [k for k, _ in items]
    
    loc_with_base = [
        (position, str(genome[position.rsplit('_', 1)[0]][int(position.rsplit('_', 1)[1])-1]).upper())
        for position in loc
    ]

    test_counts = np.array([ ## shape = (num_positions, 4) in order of VALID_BASES
        [v.get(base, 0) for base in VALID_BASES]
        for _, v in items
    ])

    counts_by_ref = defaultdict(list)
    counts_by_ref = {       ## Collect all positions into groups by referce base, then calcualte proportions of each base
        base: np.array([
            row / row.sum()  
            for (_, b), row in zip(loc_with_base, test_counts) if b == base
        ])
        for base in VALID_BASES
    }

    mean_by_ref = {     ## Calculate mean for each base in array by reference base
        base: counts_by_ref[base].mean(axis=0)
        for base in counts_by_ref
    }

    return mean_by_ref


def _shrink_position_alphas(VALID_BASES, FILE_PATHS, global_alphas, init_alphas, median_depth, genome):
    """
    Shrinks position specific alphas towards the global alpha for each reference base.
    Shrinkage is based on emperical calcualtion, where:
    Alpha_fin = (Alpha_pos * Gamma) + (Alpha_gbl * (1-Gamma))
    Gamma = Depth_pos / (Depth_pos + Tau)
    Tau is set at 0.04 * median depth of the reference. This value was found emperically. 
    """  
    first_ref = _load_reference_file(FILE_PATHS[0])
    loci = list(first_ref.keys())

    ## Tau is how much weight we give the global alphas compared to position
    ## Currently set to 4% of median value
    tau = median_depth * 0.04    

    final_alphas = {}
    for l_idx, locus in enumerate(loci):

        chromosome, pos_str = locus.rsplit('_',1)
        ref_base = str(genome[chromosome][int(pos_str)-1]).upper()

        alpha_pos = init_alphas[l_idx]
        depth = alpha_pos.sum()
        gamma = depth / (depth + tau)   ## The lower the position depth, the more we rely on the global, and vice versa
        
        alpha_gbl = np.array(global_alphas[ref_base])
        if alpha_gbl.size == 0:
            alpha_fin = alpha_pos
        else:
            alpha_fin = (gamma * alpha_pos) + ((1-gamma) * alpha_gbl)

        final_alphas[locus] = alpha_fin #* 2

    return final_alphas



def define_alphas(test_refset, genome, snp_controls):
    """
    Runs the entire script:
    - Loads files
    - Gets reference alphas for each position (counts + 1e-12), for each reference file
    - Weights reference files by depth per position; combines alphas from each reference by weighting
    - Defines global alpha-values for each reference base [ACGT]
    - Gets mean error rates of the test file for each reference base [ACGT]
    - Shifts global alphas such that the mean matches the test file
    - Shrinks position-speicifc alphas towards the mean based on position specific depth. 
    """
    with open(snp_controls, "r") as f:
        control_files = json.load(f)
    file_paths = [value for value in control_files.values()]

    ## Calculate Glboal Alphas for each base postion
    gbl_alpha_ref, position_alphas, median_depth_ref = _estimate_alphas(VALID_BASES, file_paths, genome)
    
    ## Calcualte the mean error rate of the test sample per reference base (mean by reference)    
    mbr = _calculate_test_errs(test_refset, genome)

    # Adjust global Dirichlet alphas so their mean matches the test batch,
    # while keeping total concentration (alpha0 / dispersion) unchanged.
    adjusted_global_alphas = {
        ref_base: (
            np.sum(global_alpha) * mbr[ref_base]
        )
        for ref_base, global_alpha in gbl_alpha_ref.items()
    }

    ## Shrink alphas towards the mean based on emperical shrinkage data
    final_alphas = _shrink_position_alphas(
        VALID_BASES,
        file_paths, 
        adjusted_global_alphas,
        position_alphas, 
        median_depth_ref, 
        genome
    )

    return final_alphas
