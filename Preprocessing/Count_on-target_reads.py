import pysam
import argparse


def count_reads_in_bed(bam_file, bed_file):
	"""
	Count the number of reads in a BAM file that overlap regions specified in a BED file.

	Args:
		bam_file (str): Path to the BAM file.
		bed_file (str): Path to the BED file.

	Returns:
		dict: Dictionary with regions from the BED file as keys and counts as values.
	"""
	# Open the BAM file
	bam = pysam.AlignmentFile(bam_file, "rb")
	
	# Initialize a dictionary to store counts
	region_counts = {}
	total_count = 0
	
	# Read the BED file
	with open(bed_file, "r") as bed:
		for line in bed:
			if line.startswith("#") or line.strip() == "":
				continue  # Skip comments and empty lines
			
			# Parse BED line (chrom, start, end)
			fields = line.strip().split()
			chrom = fields[0]
			start = int(fields[1])
			end = int(fields[2])
			
			# Initialize count for the region
			region_key = f"{chrom}:{start}-{end}"
			region_counts[region_key] = 0
			
			# Fetch reads overlapping the region
			for read in bam.fetch(chrom, start, end):
				region_counts[region_key] += 1
				total_count += 1
	
	# Close the BAM file
	bam.close()
	
	return region_counts, total_count


def main():
	parser = argparse.ArgumentParser(description="Count reads mapped to regions in a BED file.")
	parser.add_argument("-b", "--bam", required=True, help="Path to the BAM file.")
	parser.add_argument("-l", "--bed", required=True, help="Path to the BED file.")
	args = parser.parse_args()
	
	# Count reads in BED regions
	region_counts, total_count = count_reads_in_bed(args.bam, args.bed)
	
	# Print results
	print("Region\tCount")
	for region, count in region_counts.items():
		print(f"{region}\t{count}")
	print(f'\n{total_count}')


if __name__ == "__main__":
	main()
