---
title: Preprocessing Instructions

---

# Preprocssing Intsructions
## Extracting UMI information and pre-processing for primer trimming

Primer sequences must be removed before analysis to avoid introducing non-biological sequences that interfere with variant calling. However, this also removes the UMI, which we need. To preserve it, the UMI is extracted and added to the read ID.

The `Preprocessing.py` script performs this step and also removes duplicate reads (same ID). It requires only two positional arguments: input and outout directory locations.

```
input    The directory location containing input Fasta files
output   The directory location where files will end up. Defaults to the input directory
```

Example Usage:
```
python3 Preprocessing.py /path/to/input/ /path/to/output/
```

The output files will be in FASTA format with the suffix `_UMI.fasta`.

### Notes on UMI identification
Currently, this script can only identify UMIs from Nanopore sequenced RCA reads made from the QiaSeq Pro panel.

After BB-8 post-processing, all reads are in 5′ → 3′ orientation. Each read begins with the sample index, followed by the I5 Illumina adapter, the UMI spacer sequence, and then the UMI itself. Because of this fixed order, the UMI can be extracted from a precise location in the read. If a different library preparation technique is used, however, the wrong sequence may be extracted.

## Filtering reads by average quality and repeat number

During BB-8 consensus calling, reads have several bits of information added to the read header. The new header after post-processing will have all the information separated by underscores. We then added the UMI sequence to the end of the read ID in the previous step. In order:

| Field                 | Description  |
|-----------------------|--------------|
| Read ID               | The original read identifier |
| AverageQuality        | Average quality of the read after consensus calling |
| TotalReadLength       | Total length of the original read before consensus calling |
| NumberOfRepeats       | The total number of repeats in the entire read |
| SubReadLength         | Length of the subread, inclusive of the splint and adapters |
| Trimmed_SubReadLength | Length of the subread after post-processing adapter trimming and splint removal |
| Sequenced Strand | FORWARD or REVERSE, refers to the strand that went through the ONT pore |
| UMI | UMI-sequence extracted in the previous step |

Filtering for quality and repeat number is done using these field codes.


### Filtering by average quality
This is performed by the `Subset_by_quality_score.py` script from within the Preprocessing directory of RCA_Var. 

This script has two positional arguments, which must be the first two arguments:
```
input    The directory location containing input Fasta files
output   The directory location where files will end up. Defaults to the input directory
```

There are also several optional arguments:
```
  --help, -h            show this help message and exit
  --min_quality , -q    Minimum average Phred quality for filtered reads. Defaults to 20
  --split               Turns on split mode. This will split reads into multiple files based upon quality score
  --split_min Integer.  Sets the lower bound for splitting. Split mode only. Defaults to 20
  --split_max Integer.  Sets the upper bound for splitting. Split mode only. Defaults to 40
  --round ROUND         Rounds average quality to this numeric value e.g. `2` will round to the nearest number divisible by 2. Only useful in split mode. Defaults to 2
```

Example usage:
```
python3 Subset_by_quality_score.py /path/to/input/ /path/to/output/ -q 10
python3 Subset_by_quality_score.py /path/to/input/ /path/to/output/ --split --split_min 20 --split_max 40 --round 2
```

### Filtering by repeat number
This is performed by the `Subset_by_repeat.py` script from within the Preprocessing directory of RCA_Var. 

This script has two positional arguments, which must be the first two arguments:
```
input    The directory location containing input Fasta files
output   The directory location where files will end up. Defaults to the input directory
```

There are also several optional arguments:
```
  --help, -h            show this help message and exit
  --min_repeat , -r    Minimum average repeat count for filtered reads. Defaults to 3
  --split               Turns on split mode. This will split reads into multiple files based upon repeat number
  --split_min Integer.  Sets the lower bound for splitting. Split mode only. Defaults to 0
  --split_max Integer.  Sets the upper bound for splitting. Split mode only. Defaults to 10
```

Example usage:
```
python3 Subset_by_repeat.py /path/to/input/ /path/to/output/ -r 3
python3 Subset_by_repeat.py /path/to/input/ /path/to/output/ --split --split_min 0 --split_max 10
```

## Trimming Adapters and primer sequences
### Read Realignment and Trimming

During BB-8 post-processing, all reads were realigned in the 5' → 3' direction, starting with the Illumina I5 adapter and UMI sequence. The splint and initial adapter sequences were also trimmed at this stage; however, parts of these sequences are maintained to enable de-multiplexing and UMI identification. The Sequence in **square brackets** was used for trimming in BB-8
 
```
{Splint}AATGATA[CGGCGACCACCGAGATCTA]CAC{i5 Index} \
CTTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG{UMI+Spacer} \
{Target Read}{RSP Spacer}ATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
{i7 Index}AT[CTCGTATGCCGTCTTCTGCTTG]{Splint}
```

Therefore, we want to trim on the UMI spacer and the RSP spacer:
```
UMI Spacer (5'->3'): CATTCGAGTCAT
RSP Spacer (5'->3'): CAAAACGCAATACTGTACTGGA
```

We can run the following Cutadapt command:
```
./cutadapt -a GCATTCGAGTCAT...CAAAACGCAATACTGTACTGGAGATCGGA --untrimmed-output <UntrimmedReads.fasta> -o <OutputTrimmedFile.fasta> <InputFile.fasta>
```

By using `...` between our adapter sequences we are 'linking' them together. Reads will only be trimmed if both sequence can be found in the correct orientation. Only the internal sequence will be kept.

## Running Cutadapt to remove primer sequences
The primer binding sequence must be removed to prevent biased downstream analyses and ensure accurate representation of the target DNA. 

The `Cutadapt_primer_trimming.py` will envoke cutadapt to trim the primers from multiplexed sequencing data; only reads with a matching primer will be trimmed and all others will be added to an 'Untrimmed' file. It has two positional arguments:

```
input    The directory location containing input Fasta files to be trimmed
primers  The path to the primers.tsv file
```

There are also two optional arguments:
```
--match_len, -m    Minimum sequence length for a primer to be trimmed. Defaults to 10
--output, -o       The directory location where files will end up. Defaults to the input directory
```

Example Usage:
```
python3 Cutadapt_primer_trimming.py /path/to/input/directory/ /path/to/primers.tsv -o /path/to/output/directory -m 10
```

## Mapping
Mapping of consensus called reads is performed with minimap2. Documentation and build instructions can be found [here](https://lh3.github.io/minimap2/minimap2.html) and [here](https://github.com/lh3/minimap2/releases). Samtools is also required.

### Indexing the genome
Indexing the genome is performed with the `-d` option:

```
./minimap2 -d <genome_index.mmi> <genome.fasta>
```
There are several additional indexing options; however, is recommended to just run the default. Visit the [minimap2 documentation](https://lh3.github.io/minimap2/minimap2.html) for details.


### Mapping reads

`MapFiles.sh` will map, filter, and index multiple files: 

To run this file:
```
chmod +x MapFiles.sh
./MapFiles.sh <input dir> <output dir> <genome.fa>
```

To easily keep write out the command we can save locations of files (edit these if locations are different):
```
Input=/path/to/input/files/
Output=/path/to/output/files/
Genome=/path/to/genome.fa
Index=/path/to/index.mmi
```
Then we can just run:
`./MapFiles.sh $Input $Output $Genome $Index`

Specifically, this script:
1. Loads multiple files into a queue to be mapped
2. Maps target files and sorts the resulting BAM files with: 
    `./minimap2 -ax sr --secondary=no --cs "$Index" "$file" | samtools sort -o "$output_bam"` 
3. Adds the MD tag and removes low quality mapping with:
    `samtools calmd -u "$output_bam" "$Genome" |
        samtools view -q 20 -b -o "$output_bam_calmd"`
4. Indexes the genome with:
    `samtools index "$output_bam_calmd"`

**`MapFiles.sh` must be run from the directory containing minimap2**

## Filtering out high error reads and finalising UMIs

### Removing high error reads
Despite consensus calling, some reads (roughly 3-4%) still have a very high error rate. This appears to happen with low repeat count reads. I suspect it is the result of two or more subreads with a large amount of errors.

These reads need to be filtered out because they add unnecessary noise, and any variant found within them is questionable. This is performed with the `ctRCA_8_Removing_high_NM_reads.py` script. 

This script calculates the total number of changes (note: this is not edit distance, as insertions or deletions greater than 1 base-pair are still considered a single change). 

To run:
```
Python3 ctRCA_8_Removing_high_NM_reads.py <input-file-folder> <output-folder
```

This will take all `.bam` files in the input directory and process them, and output them in the output directory. 

### Finalising UMIs.
UMIs are initially extracted in the `ctRCA_2_Preprocessing.py` script. However, errors in the UMI sequences introduced during sequencing may split families. Furthermore, reads may have the same UMI by chance, but map to different locations. 

The `ctRCA_9_UMI_Grouping.py` script groups reads by start location (i.e. primer binding sites of the same amplicon). Reads with the same UMI, or with UMI sequences 1 hamming distance apart, are grouped together as a form of error correction.

Final UMI sequences, strand info, and primer binding sites are added to the BAM alignment as a UB:Z tag 

To run:
```
Python3 ctRCA_9_UMI_Grouping.py <input-bam> <output-bam>
```
This file only processes one at a time. 