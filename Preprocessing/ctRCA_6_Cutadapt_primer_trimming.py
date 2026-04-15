import subprocess
import os
import argparse


def parse_args():
	"""
	Parses command line arguments.
	"""
	parser = argparse.ArgumentParser(description='Wrapper for cutadapt to trim primer seqeunces from ctRCA data.',
	                                 add_help=True,
	                                 prefix_chars='-')
	# Input and Output file locations
	parser.add_argument('input', type=str, action='store',
	                    help='The directory location containing input Fasta files')
	parser.add_argument('primers', type=str, action='store',
	                    help='List of primer sequences to be trimmed. The file should only contain primer sequences each separated by a new line')
	parser.add_argument('--output', '-o', type=str, action='store', default='',
	                    help='The directory location where files will end up. Defaults to the input directory')
	parser.add_argument('--match_len', '-m', type=int, action='store', default=10,
	                    help='Minimum sequence length for a primer to be trimmed. Defaults to 10')
	
	return parser.parse_args()


def reverse_complement(sequence):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	reverse_complement_sequence = "".join(complement.get(base, base) for base in reversed(sequence))
	return reverse_complement_sequence


args = parse_args()

input_directory = args.input
if args.output == '':
	output_directory = args.input
else:
	output_directory = args.output
os.makedirs(args.output, exist_ok=True)
primer_file = args.primers

file_list = [file for file in os.listdir(input_directory) if file.endswith('.fasta') and not file.startswith('Untrimmed')]

for file in file_list:
	output_fastq = output_directory + file.rsplit('.', 1)[0] + 'PTrim.fasta'
	untrimmed_fastq = output_directory + 'Untrimmed' + file.rsplit('.', 1)[0] + 'PTrim.fasta'
	with open(primer_file) as p:
		primer_sequences = [reverse_complement(line.strip('\n').split('\t')[1].upper()) for line in p]

	command = [
		"cutadapt",
		"-O", f"{args.match_len}",
		"-o", output_fastq,
		"--untrimmed-output", untrimmed_fastq,
	]
	for primer in primer_sequences:
		command.extend(["-a", primer + 'X'])
	
	command.append(input_directory + file)
	
	try:
		result = subprocess.run(command, text=True, capture_output=True, check=True)
		print("Cutadapt ran successfully!")
		print(result.stdout)
	except subprocess.CalledProcessError as e:
		print("Error running Cutadapt:")
		print(e.stderr)
