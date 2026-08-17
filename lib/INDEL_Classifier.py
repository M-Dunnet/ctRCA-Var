import json
import ast
import numpy as np
import pandas as pd
import scipy.stats
from collections import defaultdict

class INDEL_Analyzer:
    """
    Main class for analysing INDEL variants. This class holds all the logic for:
        1. Building Beta-Binomial parameters from control datasets.
        2. Monte Carlo simulation of INDEL events to determine significance.
        3. Log-likelihood ratios for selecting alternate allele.
    """

    def __init__(
        self, 
        test_file, 
        target_indels, ## 
        control_data, 
        nn=500000, 
        eps=1e-11,
        min_depth=500,
        min_alt_count=5,
        min_alt_freq=0.001,
        pval_threshold=1.0
    ):
        ## Configuration & Parameters
        self.target_deletions_path = target_indels
        self.control_data_path = control_data
        self.test_file = test_file
        self.nn = nn
        self.eps = eps
        
        ## Filtering thresholds
        self.min_depth = min_depth
        self.min_alt_count = min_alt_count
        self.min_alt_freq = min_alt_freq
        self.pval_threshold = pval_threshold
        
        ## Data containers populated during init
        self.target_positions = self._load_target_positions()
        
        ## Load control data and cache both proportions (for beta) and absolute counts (for LLR)
        self.control_proportions, self.aggregated_ctrl_alleles = self._parse_control_data()
        
        ## Build Beta parameters
        self.beta_params = self._build_beta_parameters()

        ## Calculate PValues
        self.variant_calls = self.calculate_indel_pvalue()

    @staticmethod
    def _collapse_strand_counts(count_data):
        """
        Collapses strand counts from basecall data.
        Data is loaded in as: {'(Basecall_1, '-')': count_a, ...}
        Transforms to: {'Basecall_1': count_a + count_b}
        """
        collapsed = {}
        for position, counts_dict in count_data.items():
            collapsed[position] = {}
            for key, counts in counts_dict.items():
                base, _ = ast.literal_eval(key)
                collapsed[position][base] = collapsed[position].get(base, 0) + counts
        return collapsed

    def _load_target_positions(self):
        """Loads a CSV containing indels of interest and maps target positions."""
        dat = pd.read_csv(self.target_deletions_path)
        variant_dict = defaultdict(list)

        for pos, wt, mut in zip(
            dat["Chromosome"].astype(str) + "_" + dat["Start.Position"].astype(str),
            dat["Ref.Allele"].astype(str),
            dat["Alt.Allele"].astype(str)
        ):
            variant_dict[pos].append((wt, mut))
            
        return dict(variant_dict)

    def _parse_control_data(self):
        """
        Loads control JSON files exactly once. 
        Returns:
            proportions_by_position: List of tuples (sample_name, proportions) for Method of Moments.
            aggregated_ctrl_alleles: Dict mapping positions to their absolute control allele counts (for LLR).
        """
        with open(self.control_data_path) as f:
            ctrl_files = json.load(f)

        proportions_by_position = []
        aggregated_ctrl_alleles = defaultdict(lambda: defaultdict(int))

        for sample_name, filepath in ctrl_files.items():
            with open(filepath) as f:
                data = json.load(f)
            
            # 1. Aggregate proportions for Beta-Binomial
            collapsed = self._collapse_strand_counts(data)
            proportions = {}
            for pos, base_counts in collapsed.items():
                depth = sum(base_counts.values())
                total_dels = sum(count for base, count in base_counts.items() if base != ".")
                p = (total_dels + self.eps) / (depth + 2 * self.eps)
                proportions[pos] = (p, depth)
                
                # 2. Accumulate raw allele counts for LLR later
                for base, count in base_counts.items():
                    aggregated_ctrl_alleles[pos][base] += count

            proportions_by_position.append((sample_name, proportions))
            
        return proportions_by_position, aggregated_ctrl_alleles

    def _build_beta_parameters(self):
        """Estimates Beta-binomial distribution parameters using method of moments."""
        pos_values = defaultdict(list)
        for _, pdict in self.control_proportions:
            for pos, p in pdict.items():
                pos_values[pos].append(p[0])
        
        stats = {}
        for pos, values in pos_values.items():
            # Only process positions we actually care about
            if pos not in self.target_positions:
                continue
                
            arr = np.array(values, dtype=float)
            n = len(arr)
            if n == 0:
                continue

            mu = np.mean(arr)
            var = np.var(arr, ddof=1) if n > 1 else 0.0

            if var == 0:
                stats[pos] = {"n": n, "mean": mu, "var": 0.0, "phi": np.inf, "alpha": np.inf, "beta": np.inf}
                continue

            phi = min((mu * (1 - mu) / var) - 1, 20000)
            stats[pos] = {
                "n": n,
                "mean": mu,
                "var": var,
                "phi": phi,
                "alpha": mu * phi,
                "beta": (1 - mu) * phi
            }

        return pd.DataFrame(stats).T

    def _run_beta_binomial(self, position, alt_idx, alt_depth):
        """Generates an empirical p-value using a Monte Carlo beta-binomial simulation."""
        alpha_value = self.beta_params.loc[position]['alpha']
        beta_value = self.beta_params.loc[position]['beta']

        ran_beta = np.random.beta(alpha_value, beta_value, self.nn)
        ran_idx = np.random.binomial(alt_depth, ran_beta)

        k = np.sum(ran_idx >= alt_idx)
        pval = (k + 1) / (self.nn + 1)
        return pval

    def _select_alt_allele(self, position, obs_data):
        """Selects the most likely alternate allele based on pre-loaded control sample events."""
        test_depth = sum(obs_data.values())
        
        # Fetch pre-aggregated control counts for this specific position
        control_alleles = self.aggregated_ctrl_alleles.get(position, {})
        control_depth = sum(control_alleles.values())

        p_values = {}
        for base, obs_count in obs_data.items():
            if base == '.':
                continue

            ctrl_count = control_alleles.get(base, 0)
            expected_rate = (ctrl_count + 1) / (control_depth + 1)
            expected_count = expected_rate * test_depth

            llr = obs_count * np.log((obs_count + 1) / (expected_count + 1)) - max(obs_count - expected_count, 0)
            pval = scipy.stats.chi2.sf(llr * 2, df=1)
            p_values[base] = pval

        if not p_values:
            return None, 1.0, {}
            
        top_base = min(p_values, key=p_values.get)
        return top_base, p_values[top_base], p_values

    def calculate_indel_pvalue(self):
        """Primary function for determining significant INDEL results."""
        if isinstance(self.test_file, str):
            with open(self.test_file) as f:
                test_data = json.load(f)
        elif isinstance(self.test_file, dict):
            test_data = self.test_file
        else:
            raise ValueError(f"Deletion data is in the wrong format: {type(self.test_file)}. Expected dict or JSON filepath.")

        working_data = self._collapse_strand_counts(test_data)
        
        # Accumulate results in a list first (much faster than vcf.loc[len(vcf)])
        results = []

        for position, counts in working_data.items():
            if position not in self.beta_params.index:
                continue

            alt_depth = sum(counts.values())
            # Bug fix: use .get('.', 0) to prevent KeyError if no reference alleles exist
            alt_idx = alt_depth - counts.get('.', 0)

            # Apply configurable thresholds
            if alt_depth < self.min_depth: continue
            if alt_idx < self.min_alt_count: continue
            if alt_idx / alt_depth < self.min_alt_freq: continue

            pval = self._run_beta_binomial(position, alt_idx, alt_depth)

            if pval >= self.pval_threshold:
                continue

            allele, llr_allele, all_llr = self._select_alt_allele(position, counts)

            allele_selection_info = f"ASEL:{llr_allele}"
            allele_summary_info = f"ASUM:{','.join(f'{k}:{v}' for k, v in all_llr.items())}"
            info_column = ";".join([allele_selection_info, allele_summary_info])

            chrom, pos = position.rsplit('_', 1)
            results.append({
                'Chromosome': chrom,
                'Position': pos,
                'Type': 'INDEL',
                'Ref Allele': allele,
                'Alt Allele': '-',
                'Alt Count': alt_idx,   ## TOTAL INDEL load. 
                'Depth': alt_depth,
                'VAF': alt_idx / alt_depth,
                'P-value': pval,
                'Info': info_column
            })

        return pd.DataFrame(results)
