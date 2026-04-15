"""
This script is the last step in preprocessnig data before variant calling. 
It takes the preliminary UMI seqeunce identified in step 4: Preprocessing, 
groups reads by strand and start position (i.e. all came from the same amplicon target),
clusters UMIs within each group by Hamming distance <= 1 to account for sequencing errors,
selects a representative UMI for each cluster based on frequency and distance to other UMIs in the cluster,
and rewrites the BAM file with the representative UMI tags as UB:Z. 
"""

## Built in modules
from collections import Counter, defaultdict
import argparse
import os
import itertools
import warnings


## External modules
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", message="networkx backend defined more than once: nx-loopback")
    import networkx as nx
import pysam


def hamming_distance(seq1, seq2):
    """Calculate the Hamming distance"""
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of the same length")
    return sum(el1 != el2 for el1, el2 in zip(seq1, seq2))


def select_representative_umi(g, component):
    """
    Select a representative UMI from a cluster based on frequency, lowest total 
    Hamming distance to others or lexicographical order
    """
    ## Get read IDs and UMIs in the cluster
    cluster_reads = list(component)
    cluster_umis = [g.nodes[r]['umi'] for r in component]
    
    ## Count UMI frequencies
    umi_counts = Counter(cluster_umis)
    max_count = max(umi_counts.values())
    candidates = [umi for umi, count in umi_counts.items() if count == max_count]
    
    ## If one UMI is clearly dominant chose that one
    if len(candidates) == 1: 
        return candidates[0]
    
    ## Check that all UMIs have same frequency and there are 3 or fewer candidates, is so choose first UMI
    if len(set(umi_counts.values())) == 1 and len(candidates) <= 3: 
        return candidates[0]
    
    ## Otherwise calculate total Hamming distances and choose UMI with lowest distance to all others
    min_distance = float('inf')
    representative_umi = None
    for umi in candidates:
        total_distance = sum(hamming_distance(umi, other) for other in cluster_umis if other != umi)
        if total_distance < min_distance:
            min_distance = total_distance
            representative_umi = umi   

    return representative_umi


def generate_hamming_neighbors(umi):
    """
    Return all UMIs that differ by 1 nucleotide from the input UMI.
    """
    neighbors = []
    bases = "ACGT"
    umi_list = list(umi)
    for i, b in enumerate(umi_list):
        for new_base in bases:
            if new_base != b:
                neighbor = umi_list[:]
                neighbor[i] = new_base
                neighbors.append("".join(neighbor))
    return neighbors


def cluster_umi_groups(group_id, read_umi_pairs):
    """
    Cluster UMIs within a group by Hamming distance <= 1 to account for sequencing errors.
    """
    ## Add nodes and edges based on UMI similarity // Create a dict to map UMIs to read IDs
    g = nx.Graph()
    umi_to_reads = defaultdict(list)
    for read_id, umi in read_umi_pairs:
        g.add_node(read_id, umi=umi)
        umi_to_reads[umi].append(read_id)
    
    ## Connect all reads with identical UMIs 
    for umi, reads in umi_to_reads.items():
        for r1, r2 in itertools.combinations(reads, 2):
            g.add_edge(r1, r2)
        
        # Connect reads with UMIs differing by 1. Neighbours are generated and checked against existing UMIs
        for neighbor in generate_hamming_neighbors(umi):
            if neighbor in umi_to_reads:
                for r1 in reads:
                    for r2 in umi_to_reads[neighbor]:
                        g.add_edge(r1, r2)
    
    ## For each cluster/connected set of nodes, find the UMI that is representative of the group
    clusters = {}
    for component in nx.connected_components(g):
        cluster_reads = list(component)
        representative_umi = select_representative_umi(g, component)
        representative_umi = f"{representative_umi}_{group_id[0]}_{group_id[1]}"
        clusters[representative_umi] = cluster_reads

    return clusters


def inital_umi_collection(bam):
    """Counts UMIs in given BAM file within regions specified in BED file"""
    
    ## Init and define variables
    read_groups = defaultdict(list)
    buffer = 10  ## Allowable misalignment buffer for primer binding site
    samfile = pysam.AlignmentFile(bam, "rb")
    seen_reads = set()
                
    ## Fetch reads in region, skip unmapped reads
    for read in samfile.fetch():
        ## Skip unmapped, secondary and supplementary reads
        if read.is_unmapped:
            continue
        if read.is_secondary:
            continue
        if read.is_supplementary:
            continue

        ## Define strand, UMI and read start position
        strand = '-' if read.is_reverse else '+'
        chromosome = samfile.get_reference_name(read.reference_id)
        start = read.reference_start if read.is_reverse else read.reference_end
        umi = read.query_name.split('_')[-1]
        read_id = read.query_name

        if read_id in seen_reads:
            continue
        seen_reads.add(read_id)

        ## Group reads by strand and start position (with buffer for misalignments)
        found_key = None
        for pos in range(start - buffer, start + buffer + 1):
            key = (strand, f"{chromosome}:{start}")
            if key in read_groups:
                found_key = key
                break
        if found_key:
            read_groups[found_key].append((read_id, umi))
        else:
            read_groups[(strand, f"{chromosome}:{pos}")].append((read_id, umi))

    return read_groups


def process_bam(input_bam, output_bam):
    """
    Tag all reads with error corrected UMIs and writes to an output BAM.
    First, Gets initial UMI result -> these are grouped by strand and start position (i.e. all came from the same amplicon target).
    Read ID and UMI are stored as a tuple in the dict values -> cluster UMIs and account for errors in sequening with hamming distance
    Get a representative UMI for each cluster -> rewrite BAM file with representative UMI tags as UB:Z
    """
    ## Group reads by strand and start position and identify inital UMI seqeunces
    result = inital_umi_collection(input_bam) 

    ## Cluster UMIs in each group by hamming distance representative UMI for each family
    global_clusters = {}
    for key in result:
        clusters = cluster_umi_groups(key, result[key])
        global_clusters.update(clusters)

    ## Rewrite BAM file with representative UMI tags as UB:Z
    with pysam.AlignmentFile(input_bam, "rb") as infile, \
         pysam.AlignmentFile(output_bam, "wb", header=infile.header) as outfile:
        
        ## Invert clusters dict: read_id -> representative UMI
        read_to_umi = {read_id: umi for umi, read_ids in global_clusters.items() for read_id in read_ids}

        ## Rewrite BAM file with representative UMI tags as UB:Z
        for read in infile: 
            if read.query_name in read_to_umi:
                umi = read_to_umi[read.query_name]
                read.set_tag("UB", umi, value_type="Z")
            else:
                continue
            outfile.write(read)


def argument_parser():
    parser = argparse.ArgumentParser(description="UMI identification script",
                                     add_help=True)
    parser.add_argument("input_bam", type=str, action='store', help="Input BAM file")
    parser.add_argument("output_bam", type=str, action='store', help="Output BAM file with UMI tags added")

    return parser.parse_args()


if __name__ == "__main__":
    args = argument_parser()

    ## Check input bam
    if not args.input_bam.endswith(".bam"):
        raise ValueError("Input file must be a BAM file with .bam extension")
    
    ## Check output bam
    if not args.output_bam.endswith(".bam"):
        raise ValueError("Output file must be a BAM file with .bam extension")
    
    # output_dir = os.path.dirname(args.output_bam)
    # if not os.path.exists(output_dir):
    #     os.makedirs(output_dir, exist_ok=True)

    process_bam(args.input_bam, args.output_bam)


