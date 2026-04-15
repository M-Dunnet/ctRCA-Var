import json
import csv
import pandas as pd
import ast


def config_reader(config_in="files/config.json"):
    """
    Reads in config file. This file must be named `config.json` and must be placed within the ctRCA-Var/files/ directory.
    Alternatively, change the file path above.
    """

    with open(config_in, 'r', encoding='utf-8') as file:
        config = json.load(file)

    return config


def merge_config_with_args(config, args):
    """Merge nested JSON config with flat CLI args; CLI overrides config."""
    final_config = config.copy()

    ## Flatten CLI args for merge
    cli = {k: v for k, v in vars(args).items() if v is not None}

    ## Merge paths
    for key in final_config.get('paths', {}):
        if key in cli:
            final_config['paths'][key] = cli[key]

    ## Merge parameters
    for key in final_config.get('parameters', {}):
        if key in cli:
            final_config['parameters'][key] = cli[key]

    ## Merge plots
    for key in final_config.get('plots', {}):
        if key in cli:
            final_config['plots'][key] = cli[key]

    return final_config


def split_refset_dict(reference_dict, csv_file): ## TODO This should probably be an import?
    """
    Split a reference dictionary into two dictionaries based on positions
    present in a CSV file.

    The reference dictionary is expected to have keys in the format:
    "<chrom>_<pos>" (e.g., "NC_000017.11_58415591"), and values as nested
    dictionaries of base counts. Built from ReferenceSet module.

    The CSV file must contain 'chrom' and 'pos' columns. Positions are matched
    by separating each reference key into contig and position, and comparing
    against the (chrom, pos) pairs from the CSV.
    """
    ## Build a set of (chrom, pos) from CSV
    csv_positions = set()
    
    with open(csv_file, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            chrom = row['chrom']
            pos = row['pos']
            csv_positions.add((chrom, pos))

    ## Split the dictionary
    in_csv = {}
    not_in_csv = {}

    for key, value in reference_dict.items():
        chrom, pos = key.rsplit("_", 1)  ## safer than split("_")
        
        if (chrom, pos) in csv_positions:
            in_csv[key] = value
        else:
            not_in_csv[key] = value

    return in_csv, not_in_csv


def collapse_strand_counts(counts): ## TODO This should probably be an import?
    """
    Collapse strand-specific base counts into base-only counts.

    Converts keys like "('T', '+')" and "('T', '-')" → "T"
    and sums their counts.
    """
    collapsed = {}

    for key, value in counts.items():
        base, _ = ast.literal_eval(key)  ## safely parse "('T', '+')"
        collapsed[base] = collapsed.get(base, 0) + value

    return collapsed


def counts_to_proportions(counts):  ## TODO This should probably be an import?
    """
    Convert a dictionary of base counts into proportions.

    Example:
    {"T": 1430, "C": 1} → {"T": 0.9993, "C": 0.0007}
    """
    total = sum(counts.values())
    
    if total == 0:
        return {base: 0 for base in counts}  ## avoid division by zero

    return {
        base: value / total
        for base, value in counts.items()
    }


def refset_dict_to_vcf(test_refset, genome):       ## TODO This should probably be an import?
    """
    Convert collapsed reference dictionary into a long-format DataFrame.

    Parameters
    ----------
    test_refset : dict
        {position: {base: count}}
    genome : pyfaidx.Fasta or similar
        Reference genome object supporting genome[contig][pos]

    Returns
    -------
    pd.DataFrame
        Columns: Contig, Position, Ref, Alt, Depth, Alt_Count, Alt_Prop
    """
    rows = []
    for key, counts in test_refset.items():
        # Split key
        contig, pos = key.rsplit("_", 1)
        pos = int(pos)

        # Reference base (0-based indexing)
        ref_base = str(genome[contig][pos - 1]).upper()

        # Depth
        depth = sum(counts.values())
        if depth == 0:
            continue

        # Iterate over non-reference bases
        for base, count in counts.items():
            if base == ref_base:
                continue

            rows.append({
                "Contig": contig,
                "Position": pos,
                "Ref": ref_base,
                "Alt": base,
                "Depth": depth,
                "Alt_Count": count,
                "Alt_Prop": count / depth
            })

    return pd.DataFrame(rows)


def load_bed(bed_path):
    """
    Loads in BEDfile as PD.DataFrame from a path string. 
    Bed file must only contin contig, start, end, and gene columns seperated by tabs.
    """
    ## Load in BED file as PD.DataFrame
    bed_df = pd.read_csv(
        bed_path,
        sep="\t",
        comment="#",
        header=None,
        names=["Contig", "Start", "End", "Gene"]
    )

    ## Ensure positions are integers
    bed_df["Start"] = bed_df["Start"].astype(int)
    bed_df["End"] = bed_df["End"].astype(int)

    return bed_df


def annotate_gene(row, bed_df):
    """
    Takes contig and postion information for each candiate variant 
    and checks which gene that is associeted with in the BEDfile
    """
    ## Filter BED to same contig
    contig_bed = bed_df[bed_df["Contig"] == row["Contig"]]
    ## Check if position is within start-end
    match = contig_bed[
        (row["Position"] >= contig_bed["Start"]) &
        (row["Position"] <= contig_bed["End"])
    ]
    if not match.empty:
        return match["Gene"].values[0]  ## take first match
    
    return None