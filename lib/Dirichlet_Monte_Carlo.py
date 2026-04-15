import sys
import pandas as pd
import json
import numpy as np


class Variant:
	"""Represents a single variant with relevant attributes."""
	
	def __init__(self, chrom, pos, ref, alt, alt_count, depth):
		self.chrom = chrom # Chromosome
		self.pos = pos # Position in the genome
		self.ref = ref # Ref base
		self.alt = alt	# Alt base
		self.alt_count = alt_count # Count of alt bases
		self.depth = depth # Total depth of coverage (after UMI consensus calling)
		self.location = f"{chrom}_{pos}"


class DirichletSimulator:
	"""Handles Dirichlet distribution generation and multinomial sampling."""
	
	def __init__(self, nn=500000):
		self.nn = nn  # Number of Dirichlet samples
	
	def generate_dirichlet(self, alphas):
		"""Generates Dirichlet-distributed probabilities."""
		return np.random.dirichlet(list(alphas), self.nn)
	
	def generate_multinomial_samples(self, dirichlet_samples, depth):
		"""Generates multinomial samples based on Dirichlet output."""
		return [np.random.multinomial(depth, d, 1)[0] for d in dirichlet_samples]


class VariantAnalyzer:
	"""Main class for handling variant data, error data, and computing p-values."""
	
	def __init__(self, var_file, ref_file, nn=500000):
		self.var_file = var_file
		self.ref_file = ref_file
		self.nn = nn
		self.valid_bases = {"A", "C", "G", "T"}
		self.filtered_variants = self.filter_variants()
		self.simulator = DirichletSimulator(nn)
	
	def filter_variants(self):	#TODO Fix this so CONTIG_POS matches VARIANT LOCATIONS. // fix this method
		"""Filters variant locations based on available error data."""
		variant_locations = self.var_file["Contig"].astype(str) + "_" + self.var_file["Position"].astype(str)
		valid_error_keys = set(self.ref_file.keys()) & set(variant_locations)
		return self.var_file[variant_locations.isin(valid_error_keys)]
	
	def process_variant(self, variant):
		"""Processes a single variant to compute p-value."""
		
		## Get alphas from alpha dictonary
		## Alphas is a 1D array: array([A,C,G,T])
		alphas = self.ref_file[variant.location]
		alpha_mapping = {
			'A':0,
			'C':1,
			'G':2,
			'T':3
		}

		# Get the error probability of the current variant from the reference set 
		err_prob = alphas[alpha_mapping[variant.alt]] / sum(alphas)
		
		# Generate Dirichlet samples and multinomial data
		dirichlet_samples = self.simulator.generate_dirichlet(alphas)
		multinomial_samples = self.simulator.generate_multinomial_samples(dirichlet_samples, variant.depth)
		
		# Get counts for alt base from multinomial samples
		alt_idx = alpha_mapping[variant.alt]
		alt_rand = [sample[alt_idx] for sample in multinomial_samples]
		
		# Compute empirical p-value
		p_value = (sum(np.array(alt_rand) > variant.alt_count) + 1) / (self.nn + 1)
		# print(variant.chrom, variant.pos, variant.ref, variant.alt, variant.alt_count, variant.depth, p_value)
		return [variant.location, variant.ref, variant.alt, variant.alt_count, variant.depth, err_prob, p_value]
	
	def run_analysis(self):
		"""Runs the entire analysis and computes empirical p-values."""
		results = []
		for _, row in self.filtered_variants.iterrows():
			variant = Variant(row["Contig"], row["Position"], row["Ref"], row["Alt"], row["Alt_Count"], row["Depth"])
			results.append(self.process_variant(variant))
		return results



# Convert results to DataFrame for easier handling
# columns = ["Location", "Ref", "Alt", "Alt_Count", "Depth", "Err_Base", "Err_Prob", "P_Value"]
# emp_pv_df = pd.DataFrame(empirical_p_values, columns=columns)
# emp_pv_df.to_csv(f"/Users/dunmi18p/Desktop/cfDNA_beta_conf/{var_file.split('/')[-1]}_empirical_p_values.tsv", sep="\t", index=False)
