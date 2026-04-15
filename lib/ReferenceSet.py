'''
Holder for error profile Class
'''

from collections import defaultdict, Counter 
import pysam
from pyfaidx import Fasta
from abc import ABC, abstractmethod

class ReferenceSetBase(ABC):
    '''Base class for reference sets'''
    def __init__(self, file_path, bed_path, genome_path, positions_to_remove=None, depth_threshold=100, collapse_umi=True):
        self.depth_threshold = depth_threshold
        self.file = file_path
        with open(bed_path) as b:
            self.bedfile = [line.strip() for line in b if not line.startswith('#')]
        self.genome = Fasta(genome_path)
        self.collapse_umi = collapse_umi
        self.reference_set = self._generate_reference_set()
        
        try:
            if positions_to_remove is None:
                self.exclude = set()
            else:
                with open(positions_to_remove) as f:
                    self.exclude = [
                        (chrom, int(start), int(end))
                        for line in f
                        if not line.startswith('#')
                        for chrom, start, end in [line.strip().split('\t')]
                    ]

        except ValueError as e:
            raise ValueError("Error: At least one value in the positions_to_remove list could not be coerced into an integer") from e

        except FileNotFoundError as e:
            raise FileNotFoundError(f"Error: The file '{positions_to_remove}' was not found.") from e
        
    def _generate_reference_set(self):
        """Generalized method to process the BAM file."""
        samfile = pysam.AlignmentFile(self.file)
        process = self._process_pileupread
        consensus_basecall_dict = defaultdict(list)
        bed_regions = [(line.split()[0], int(line.split()[1]), int(line.split()[2])) for line in self.bedfile]
        bed_regions.sort()

        for i, (chrom, start, end) in enumerate(bed_regions, start=1):
            print(f"\tBuilding reference set... Working on region {i} of {len(bed_regions)} - {round(i/len(bed_regions) * 100, 2)}% ({chrom}:{start}-{end})", ' ' * 10, end='\r')
            for pileupcolumn in samfile.pileup(chrom, start, end, max_depth=25000, truncate=True): 

                chromosome = pileupcolumn.reference_name
                position = pileupcolumn.pos + 1  # 1-based position
                umi_dict = defaultdict(list)

                ## Process each pileup read
                if self.collapse_umi:
                    for pileupread in pileupcolumn.pileups:
                        aln = pileupread.alignment  ## Avoids repeated attribute access, reuse aln instead of repeated pileupread.alignment
                        try:
                            umi = aln.get_tag('UB')
                        except KeyError:
                            continue
                        updated_position, processed_call, strand = process(pileupread, chromosome, position)
                        if processed_call is not None:
                            umi_dict[umi].append((processed_call, strand))  # Store call with strand information

                    ## UMI consensus base call
                    umi_position_consensus = self._get_consensus_base(umi_dict, chrom, position) ## TODO need argument to set family size.
                    for umi_seq in umi_position_consensus.values():
                        if updated_position is None:
                            updated_position = position
                        consensus_basecall_dict[f"{chromosome}_{updated_position}"].append(umi_seq)

                # If not collapsing UMIs, process each read uniquely
                else:
                    for pileupread in pileupcolumn.pileups:
                        updated_position, processed_call, strand = self._process_pileupread(pileupread, chromosome, position)
                        if processed_call is not None:
                            consensus_basecall_dict[f"{chromosome}_{updated_position}"].append((processed_call, strand))  # Store call with strand information
        
        print()  # For newline after progress output
        reference_set = {key: Counter(value) for key, value in consensus_basecall_dict.items()}
        return {genome_loc: {str(basecall): count for basecall, count in num_calls.items()} 
                for genome_loc, num_calls in reference_set.items() 
                if sum(num_calls.values()) >= self.depth_threshold}
    

    @abstractmethod
    def _process_pileupread(self, pileupread, chromosome, position):
        """Extracts SNP, deletion, or insertion-specific information from a pileup read."""
        pass

    @staticmethod
    def _get_consensus_base(sequence_dict, chrom, position, reference_genome, umi_family_size=1): ## TODO Might Need to Remove This
        ref_base = str(reference_genome[chrom][position - 1]).upper()  # `position` is the 1-based index, Pysam and Pyfaidx use 0-based indexing
        
        consensus_dict = {}
        for umi, read_info in sequence_dict.items():
            calls = [c for c, _ in read_info]
            strand = [s for _, s in read_info]
            total_calls = len(calls)
            consensus_threshold = total_calls * 0.65
            base_counts = Counter(calls)
            strand_counts = Counter(strand)
            
            if total_calls < umi_family_size:
                continue  # Skip UMIs that don't meet the family size threshold
            consensus_base = [base for base, count in base_counts.items() if count >= consensus_threshold]
            consensus_strand = [s for s, count in strand_counts.items() if count >= consensus_threshold]

            if len(consensus_strand) == 0:
                continue
            if len(consensus_base) == 0:
                consensus_base = str(ref_base).split()[-1]
            else:
                consensus_base = consensus_base[0]
            consensus_dict[umi] = (consensus_base, consensus_strand[0])

        return consensus_dict
    
    def filter_positions(self):  # This will need changing if I want to preserve the original JSON file
        for location in list(self.reference_set.keys()):
            chrom, position = location.rsplit('_', 1)
            if self.exclude is not None and any(
                    exclude_chrom == pileupcolumn.reference_name and exclude_start <= pileupcolumn.pos + 1 <= exclude_end
                    for exclude_chrom, exclude_start, exclude_end in self.exclude
                ):
                    continue
    

class SNP_ReferenceSet(ReferenceSetBase):
    """Handles SNP reference set generation."""
    
    def _process_pileupread(self, pileupread, chromosome, position):
        """Extracts SNP base call from pileup read."""
        aln = pileupread.alignment
        if aln.is_supplementary or aln.is_secondary:
            return None, None, None
        if pileupread.is_refskip or pileupread.is_del or pileupread.indel != 0:
            return None, None, None
        if pileupread.query_position is None:
            return None, None, None
        
        if aln.is_forward:
            strand = '+'
        else:
            strand = '-'
        return position, str(aln.query_sequence[pileupread.query_position]).upper(), strand
    
    def _get_consensus_base(self, sequence_dict, chrom, position):
        return super()._get_consensus_base(sequence_dict, chrom, position, reference_genome=self.genome)
    

class Del_ReferenceSet(ReferenceSetBase):
    """Handles Deletion reference set generation."""
    
    def _process_pileupread(self, pileupread, chromosome, position):
        """Extracts deletion event from pileup read."""
        aln = pileupread.alignment
        if aln.is_supplementary or aln.is_secondary:
            return None, None
        if pileupread.is_refskip or pileupread.is_del:
            return None, None
        if pileupread.query_position is None:
            return None, None, None
        if aln.is_forward:
            strand = '+'
        else:
            strand = '-'
        if pileupread.indel < 0:  # Deletion event
            deletion_length = -pileupread.indel
            return position + 1, str(self.genome[chromosome][position:position + deletion_length]).upper(), strand
        return position + 1, '.', strand

    def _get_consensus_base(self, sequence_dict, chrom, position): 
        consensus_dict = {}
        for umi, read_info in sequence_dict.items():
            calls, strand = zip(*read_info)
            base_counts = Counter(calls)
            strand_counts = Counter(strand)
            total_calls = len(calls)
            consensus_threshold = total_calls * 0.65

            consensus_base = [base for base, count in base_counts.items() if count >= consensus_threshold]
            consensus_strand = [strand for strand, count in strand_counts.items() if count >= consensus_threshold]

            if len(consensus_strand) == 0:
                continue
            if not consensus_base:
                consensus_base = '.'
            else:
                consensus_base = consensus_base[0]
            consensus_dict[umi] = (consensus_base, consensus_strand[0])

        return consensus_dict


class Ins_ReferenceSet(ReferenceSetBase):
    """Handles Insertion reference set generation."""
    
    def _process_pileupread(self, pileupread, chromosome, position):
        """Extracts insertion event from pileup read."""
        aln = pileupread.alignment
        if aln.is_supplementary or aln.is_secondary:
            return None, None, None
        if pileupread.is_refskip or pileupread.is_del:
            return None, None, None
        if pileupread.query_position is None:
            return None, None, None
        if aln.is_forward:
            strand = '+'
        else:
            strand = '-'
        if pileupread.indel > 0:  # Insertion event
            insertion_length = pileupread.indel
            return position, str(aln.query_sequence[pileupread.query_position:pileupread.query_position + insertion_length + 1]).upper(), strand
        return position, '.', strand

    def _get_consensus_base(self, sequence_dict, chrom, position): 
        consensus_dict = {}
        for umi, read_info in sequence_dict.items():
            calls, strand = zip(*read_info)
            base_counts = Counter(calls)
            strand_counts = Counter(strand)
            total_calls = len(calls)
            consensus_threshold = total_calls * 0.65

            consensus_base = [base for base, count in base_counts.items() if count >= consensus_threshold]
            consensus_strand = [strand for strand, count in strand_counts.items() if count >= consensus_threshold]
            
            if len(consensus_strand) == 0:
                continue    
            if not consensus_base:
                consensus_base = '.'
            else:
                consensus_base = consensus_base[0]
            consensus_dict[umi] = (consensus_base, consensus_strand[0])

        return consensus_dict
    
