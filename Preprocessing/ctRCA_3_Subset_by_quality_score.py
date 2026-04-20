import argparse
from pyfaidx import Fasta
import os
import sys
import platform


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
	parser.add_argument('--min_quality', '-q', type=int, action='store', default=15,
	                    help='Minimum average Phred quality for filtered reads. Defaults to 15')
	parser.add_argument('--split', action='store_true', default=False,
	                    help='Turns on split mode. This will split reads into multiple files based upon quality score')
	parser.add_argument('--split_min', type=int, action='store', default=20,
	                    help='Integer. Sets the lower bound for splitting. Split mode only. Defaults to 20')
	parser.add_argument('--split_max', type=int, action='store', default=40,
	                    help='Integer. Sets the upper bound for splitting. Split mode only. Defaults to 40')
	parser.add_argument('--round', type=int, action='store', default=2,
	                    help='Rounds average quality to this numeric value e.g. `2` will round to the nearest number divisible by 2. Only useful in split mode. Defaults to 2')
	
	return parser.parse_args()


def round_to_nearest(num, x):
	"""
	Round a number to the nearest multiple of a given base.
	Parameters:
		num (int): The number to be rounded.
		x (int): The base to which 'num' will be rounded to the nearest multiple.
	Returns:
		int: The number rounded to the nearest multiple of 'x'.
	Raises:
		ValueError: If 'x' is zero (division by zero is not allowed).
	"""
	return round(num / x) * x


def quality_cutoff(args, file_path, filename):
	batch_size = 10000
	batch = []
	
	fasta_file = Fasta(file_path)
	num_reads = len(fasta_file.keys())

	with open(f"{args.output}{filename}_{str(args.min_quality)}Plus_quality.fasta", "w") as f:
		for i, (record_id, sequence) in enumerate(fasta_file.items()):
			# The Record ID contains multiple components generated after R2C2 consensus calling,
			# summarizing various read metrics. Each component is separated by an underscore (_).
			# Format:
			#   ReadID_AvgQuality_TotalLength_RepeatCount_FullSubreadLength_TrimmedSubreadLength
			
			print(f"\tWorking on {filename}: {round((i/num_reads)*100, 2)}% complete", end="\r")
			
			id_components = record_id.split('_')
			qual = float(id_components[1])
			if qual < args.min_quality:
				continue
			
			batch.append(f">{record_id}\n{sequence}\n")
			if len(batch) >= batch_size:
				f.writelines(batch)
				batch = []
				
		if batch:
			f.writelines(batch)
			

def split_files(args, file_path, filename):
	fasta_file = Fasta(file_path)
	num_reads = len(fasta_file.keys())
	
	qual_range = range(args.split_min, args.split_max + 1, args.round)
	split_qual_dict = {key: [] for key in qual_range}
	
	for i, (record_id, sequence) in enumerate(fasta_file.items()):
		# The Record ID contains multiple components generated after R2C2 consensus calling,
		# summarizing various read metrics. Each component is separated by an underscore (_).
		# Format:
		#   ReadID_AvgQuality_TotalLength_RepeatCount_FullSubreadLength_TrimmedSubreadLength
		
		print(f"\tWorking on {filename}: {round((i / num_reads) * 100, 2)}% complete", end="\r")
		
		id_components = record_id.split('_')
		qual = float(id_components[1])
		rounded_qual = round_to_nearest(qual, args.round)
		
		if rounded_qual not in split_qual_dict.keys():
			continue
		split_qual_dict[rounded_qual].append(f">{record_id}\n{sequence}\n")
	
	for key, value in split_qual_dict.items():
		with open(f"{args.output}{filename}_{key}_quality.fasta", "w") as f:
			f.writelines(value)


def main(args):
	file_list = []
	for file in os.listdir(args.input):
		if file.endswith('.fasta') or file.endswith('fa'):
			file_list.append(file)
	
	# Check file list contains FASTA files
	if len(file_list) == 0:
		print('There are no Fasta files in input folder...')
		sys.exit(1)
	
	print(f"----------------------------------------------------------\n"
	      f"Subsetting Consensus called reads by average quality score\n"
	      f"----------------------------------------------------------\n"
	      f"{len(file_list)} file(s) provided")
	
	for file in file_list:
		file_path = args.input + file
		if args.split:
			split_files(args, file_path, file)
			print('')
		else:
			quality_cutoff(args, file_path, file)
			print('')
	print('\n')


if __name__ == '__main__':
	system_platform = platform.system()
	if system_platform == 'Windows':
		os.system('cls')
	else:
		os.system('clear')
	
	## Parse arguments
	args = parse_args()
	
	## Check an input directory was provided
	if not args.input:
		print('Input folder (--input/-i) is required', file=sys.stderr)
		sys.exit(1)
	
	## Check the input folder ends with a backslash
	if not args.input.endswith('/'):
		args.input = args.input + '/'
	
	## If an output folder is not supplied, make the input folder the output
	if not args.output:
		print('An output location was not provided. Output files will be placed in the input folder')
		args.output = args.input

	## Check the output folder ends with a backslash
	if not args.output.endswith('/'):
		args.output = args.output + '/'
	
	## Create output directory if it does not exist
	os.makedirs(args.output, exist_ok=True)
	
	main(args)
