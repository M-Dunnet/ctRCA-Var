# ctRCA-Var
ctDNA variant caller for RCA consensus called data. 

This pipeline assumes that target enrichment of ctDNA moleucles has been performed with the QIASeq-pro panel, that circularisation + RCA was performed using the R2C2 protocol, and that consensus calling was peformed with either C3POa or BB-8

## Requirements
Mutliple supporting files are needed to run ctRCA-Var:
- Genome file in fasta format
- Target regions in BED format
- Exon locations of target genes in BED format
- Target gene coding sequence in fasta format
- List of known variants in tsv format
- primer binding locations in BED format
- Primer sequnces in tsv format

A minimap2 executable are requried to run the full pipeline, as well as several pip installables.

### Obtaining supporting files
See `obtaining supporting files` MarkDown for information

### Pip installables
The following modules are required:
- cutadapt
- numpy
- pandas
- pysam
- pyfaidx
- matplotlib
- seaborn

Install these requirements with pip:

```
python3 -m venv venv
source venv/bin/activate
venv/bin/pip --upgrade pip
venv/bin/pip install cutadapt
venv/bin/pip install numpy
venv/bin/pip install pandas
venv/bin/pip install pysam
venv/bin/pip install pyfaidx
venv/bin/pip install matplotlib
venv/bin/pip install seaborn
```

### minimap2 installation
To Download, visit the releases page of the GitHub repository, and download the pre-compiled binary (https://github.com/lh3/minimap2/releases).
Before minimap2 can be used, the target genome needs to be indexed:

```
./minimap2 -d <genome_index.mmi> <genome.fasta>
```

## Pipeline overview

This pipeline is designed for consensus called nanopore reads from BB-8 (or C3POa) as input. ctRCA-Var outputs mutation annotation format tsv files at the end of the pipeline.
The preprocessing file is as follows. Usage instructions for each step in this process are in the `Preprocessing Instructions` MarkDown file
```mermaid
flowchart LR
    A[Inital UMI assignment] --> B[Quality filtering]
    B --> C[Repeat filtering]
    C --> D[Adapter trimming]
    D --> E[Primer trimming]
    E --> F[Mapping]
    F --> H[High error read removal]
    H --> I[Read UMI grouping]
```
After preprocessing, BAM files are ready for ctRCA-Var variant calling

## ctRCA-Var Usage
ctRCA-Var can be run using command line arguments or through a config file. Command line arguments take precedence over the config file; however, if the same genome, BED file, etc are consistanely re-used it is helpful to add them to the config to prevent excessive commandline arguments.

Example Usage:
`python your_script.py --input <input_directory> [options]`
`--input` is required, and points to a directory (NOT a file). All files with the `bam` suffix will be processed.

ctRCA-Var has several options:
| Argument                | Value                | Description                                                                         |
| ----------------------- | -------------------- | ----------------------------------------------------------------------------------- |
| --input                 | ./bam_files          | The directory location containing input BAM files                                   |
| --output                | ./results            | The directory location where files will end up (defaults to input directory)        |
| --genome                | ./ref/genome.fa      | Path to the reference genome                                                        |
| --bedfile               | ./targets.bed        | Path to a BED file; only loci in this file will be analysed                         |
| --control_refset        | ./controls/          | Path to control sample refset files; required for Dirichlet-multinomial analysis    |
| --exon_bed              | ./ref/exons.bed      | Path to exon boundary BED file; used to determine if variants are in coding regions |
| --cds_fasta             | ./ref/cds.fa         | Path to canonical CDS sequences                                                     |
| --target_variants       | ./filters/variants.tsv | Path to target variants filter file                                                 |
| --UMI_cutoff            | 3                    | Minimum number of reads in a UMI family; families below this are ignored            |
| --UMI_proportion        | 0.75                 | Proportion of reads in a UMI group containing a mutation (0.0–1.0)                  |
| --strand_bias_threshold | 0.01                 | Significance threshold for strand bias using Fisher's exact test                    |
| --minimum_coverage      | 100                  | Minimum read coverage required when determining relevant mutations                  |
| --minimum_vaf           | 0.005                | Minimum variant allele frequency                                                    |
