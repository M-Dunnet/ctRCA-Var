---
title: Obtaining Supporting files

---

# Obtaining Supporting files
## Summary
This workflow will explain how to get the supporting files used in ctRCA analysis pipeline. The files required are:
| File                     | Used for                          | Process          | File Type |
|---------------------------|-----------------------------------|------------------|-----------|
| Target regions            | ctRCA-Var                        | Variant calling  | BED       |
| Genome sequence           | minimap2, ctRCA-Var, Cutadapt | Multiple         | FASTA     |
| Target gene exon locations     | ctRCA-Var                        | Variant calling  | BED       |
| Target gene CDS sequences      | ctRCA-Var                        | Variant calling  | FASTA     |
| Target Variants           | ctRCA-Var                        | Variant calling  | CSV       |
| Primer binding locations  | Cutadapt                         | Primer trimming  | BED       |
| Primer sequences          | Cutadapt                         | Primer trimming  | TSV       |

---

## Step 1: Target Regions BED File
This file tells RCA_Var which regions to examine. Only variants within the regions specified in the BED file will be included in the final VCF file.

The RCA_Var pipeline assumes that the QIASeq Pro Panel has been used to generate data. The regions of interest (ROI) BED file provided by QIAGEN can be used directly.

* Ensure that the genome coordinates and contig names match the genome build used for mapping.
* The final BED file should follow this structure, with fields separated by tabs.
* Lines beginning with `#` will be ignored.

The final file should be in this structure, separated by tabs:
| Chromosome    | Start      | End        | Gene |
|--------------|-----------|-----------|------|
| NC_000001.11 | 114713889 | 114713930 | NRAS |
| NC_000001.11 | 114716103 | 114716144 | NRAS |

---

## Step 2: Genome Fasta File
The genome FASTA file is used for mapping, helps RCA var determine the reference base for variants, and is used to get primer seqeunces from the primer bed file.

**The reference genome must be the same build used throughout.**

Reference genomes can be downloaded from the [NCBI website](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/). It is recommended to use the most up to date genome build. 


## Step 3: Exon locations BED file
This is used by ctRCA_Var to determine coding and non-coding variants, as well as coding sequence changes.

Data is obtained from the UCSC Genome Browser's Table Viewer. Exon coordinates will be generated for canonical transcripts only. Splice variants will not be included.

The process requires getting the canonical mRNA transcript coordinates, exon coordinates from the canonical transcripts, and joining these two datasets together.

**Getting the canonical transcript coordinates**
- Head to the UCSC genome browser website: https://genome.ucsc.edu/
- Navigate to **Tools** -> **Table Viewer**
- Under select dataset, select the desired orgnaism and genome assembly:
    ```
    clade: Mammal
    genome: Human
    assembly: GRCH38/hg38
    ```
 - Select the following dataset:
    ```
    group: Genes and Gene Predictions
    track: GENCODE V44 (or most up to date version)
    table: knownCanonical
    ```
- Select the entire genome as the region of interest:
    ```
    region: genome
    ```
- Specifcy output data format:
    ```
    output format: selected fields from primary and related tables
    output filename: UCSC_canonical.bed
    output field separator: tsv
    ```
- Click `get output`
- From `Select Fields from hg38.knownCanonical` select:
    ```
    chrom
    chromStart
    chromEnd
    transcript
    ```
- From `hg38.kgXref fields` select `geneSymbol`
- click `get output` (under `hg38.knownCanonical`)

**Getting exon coordinates**
- Head to the UCSC genome browser website: https://genome.ucsc.edu/
- Navigate to **Tools** -> **Table Viewer**
- Under select dataset, select the desired orgnaism and genome assembly. In my case this was:
    ```
    clade: Mammal
    genome: Human
    assembly: GRCH38/hg38
    ```
 - Select the following dataset:
    ```
    group: Genes and Gene Predictions
    track: GENCODE V44 (or most up to date version)
    table: knownGene
    ```
- Select the entire genome as the region of interst, unless your requirements are different
    ```
    region: genome
    ```
- Specifcy output data format:
    ```
    output format: BED
    output filename: UCSC_exons.bed
    output field separator: tsv
- Click `get output`
- Under `Create one BED record per` select `Exons` or `Coding exons` (depending on what you're after). Leave 0 bases at the end.
- Click `get bed`

**Joining Exon data and canonical transcript information**
You should now have two downloaded files:
- UCSC_exons.bed
- UCSC_canonical.bed 

Modify the exon file to seperate the transcript name from the rest of the information:
```
awk '{split ($4,a,"_"); {print $1"\t"$2"\t"$3"\t"a[1]"\t"a[3]"\t"$6}}' UCSC_exons.bed > UCSC_exons_modif.bed
```
Join the two files based on transcript identifier:
```
join -1 4 -2 4 <(sort -k4 UCSC_exons_modif.bed ) <(sort -k4 UCSC_canonical.bed) | awk '{print $2"\t"$3"\t"$4"\t"$10"\t"$1"\t"$5"\t"$6}'' | bedtools sort -i "-" > UCSC_exons_modif_canonical.bed
```

You will get a file called `UCSC_exons_modif_canonical.bed`. It will look like this:
					

| Chromosome | Start  | End   | Gene    | Transcript ID      | Exon number | Strand |
|------------|--------|-------|---------|--------------------|-------------|--------|
| chr1       | 11868  | 12227 | DDX11L2 | ENST00000000233.10 | 0           | +      |
| chr1       | 12009  | 12057 | DDX11L1 | ENST00000000233.10 | 0           | +      |
| chr1       | 12178  | 12227 | DDX11L1 | ENST00000000233.10 | 1           | +      |
| chr1       | 12612  | 12721 | DDX11L2 | ENST00000000233.10 | 1           | +      |
| chr1       | 12612  | 12697 | DDX11L1 | ENST00000000233.10 | 2           | +      |
| chr1       | 12974  | 13052 | DDX11L1 | ENST00000000233.10 | 3           | +      |
| chr1       | 13220  | 13374 | DDX11L1 | ENST00000000233.10 | 4           | +      |
| chr1       | 13220  | 14409 | DDX11L2 | ENST00000000233.10 | 2           | +      |

If you only want the exons of a single gene, use this command:
```bash
awk -F'\t' '$4 == "GENE" { print }' UCSC_exons_modif_canonical.bed > output.bed
```
Where `GENE` is the gene ID of your favorite gene.

To subset specific genes, create a case-sensitive text file (e.g., genes.txt) with one gene name per line, then run this command:

```bash
awk -F'\t' 'NR==FNR { genes[$1]; next } $4 in genes' genes.txt UCSC_exons_modif_canonical.bed > output.bed
```

## Step 4: Getting CDS FASTA sequences
This step requires three files:
* The canonical exons BED file (final file from step 3)
* The CDS sequence in FASTA format
* A table that cross-references canonical CDS sequence with HUGO gene IDs.

Once those files are ready, the transcript ID, Gene ID, and CDS ID need to be added to the CDS FASTA header for each CDS transcript. This is performed with the `Add_transcriptID_to_CDS_FASTA.py` 

**Getting the CDS FASTA Seqeunce**
- Head to the UCSC genome browser website: https://genome.ucsc.edu/
- Navigate to **Tools** -> **Table Viewer**
- Under select dataset, select the desired orgnaism and genome assembly. In 
- my case this was:
    ```
    clade: Mammal
    genome: Human
    assembly: GRCH38/hg38
    ```
 - Select the following dataset:
    ```
    group: Genes and Gene Predictions
    track: CCDS
    table: ccdsGene
- Select the entire genome as the region of interst, unless your requirements are different
    ```
    region: genome
    ```
- Specifcy output data format:
    ```
    output format: Sequence
    output filename: UCSC_Transcript_CDS_sequence.bed
    ```
- Click `get output`
- Click `submit` for `genome` - its the only option
- Ensure that only `CDS Exons`, `One FASTA record per gene`, and `All upper case` are checked, then click `get sequence`


**Getting tables to cross-reference the CDS ID with Gene ID**

- Head to the UCSC genome browser website: https://genome.ucsc.edu/
- Navigate to **Tools** -> **Table Viewer**
- Under select dataset, select the desired orgnaism and genome assembly. In my case this was:
    ```
    clade: Mammal
    genome: Human
    assembly: GRCH38/hg38
    ```
 - Select the following dataset:
    ```
    group: Genes and Gene Predictions
    track: GENCODE V44 (or most up to date version)
    table: knownGene
    ```
- Select the entire genome as the region of interst, unless your requirements are different
    ```
    region: genome
    ```
- Specifcy output data format:
    ```
    output format: selected fields from primary and related tables
    output filename: Transcript_CCDS_ID.bed
    output field separator: tsv
- Click `get output`
- Ensure the following fields are selected from `Select Fields from hg38.knownGene`:
    ```
    name
    ```
- Ensure the following fields are selected from `Linked tables`:
    ```
    ccdsInfo
    kgXref
    ```
- Scroll to the bottom of the page and click `allow selection from linked tables`. This will refresh the page, now with more table options.
- From `hg38.ccdsInfo fields` ensure `ccds` is selected. 
- from `kgXref` ensure `geneSymbol` is selected.
- These are the only two feilds required, however, these others may be useful:
    - from `Select Fields from hg38.knownGene`:
        ```
        chrom
        strand
        txStart
        txEnd
        ```
    - from `hg38.ccdsInfo fields`:
        ```
        mrnaAcc
        protAcc
        ```
    - from `hg38.kgXref fields`:
        ```
        mRNA
        spID
        refseq
        protAcc
        ```
- Additional fields will help with context, however, they will change the order of columns. Ensure you know the column number of `name` and `ccdsInfo` before continuing. 
- click `get output`
- I have chosen to have all additional columns listed above from `Select Fields from hg38.knownGene`, but none from the other fields. This gives me a table like: 

| #hg38.knownGene.name   | hg38.knownGene.chrom | hg38.knownGene.strand | hg38.knownGene.txStart | hg38.knownGene.txEnd | hg38.ccdsInfo.ccds | hg38.kgXref.geneSymbol |
|-------------------------|----------------------|------------------------|-------------------------|----------------------|---------------------|------------------------|
| ENST00000256078.10      | chr12                | -                      | 25205245                | 25250929             | CCDS8703.1          | KRAS                   |
| ENST00000311936.8       | chr12                | -                      | 25205245                | 25250929             | CCDS8702.1          | KRAS                   |
| ENST00000686877.1       | chr12                | -                      | 25205249                | 25250908             | CCDS8702.1          | KRAS                   |

**Note that ccdsInfo fields only link to the consensus cds seqeunce, knownGene fields are individual transcript variants**

**Putting it all together**
Run the `Add_transcript_to_CDS_fasta.py`. You will get an output file called `Final_CDS_sequences.fasta` ready for use. 

## Primer binding location BED file and primer sequences

Qiagen do not provide us with the exact primer sequence, rather they give us a BED file of primer binding locations and expected coverage based on read length. For example:

| Chromosome | Start      | End        | Name         | Score | Strand | ThickStart | ThickEnd   |
|------------|------------|------------|--------------|-------|--------|------------|------------|
| chr1       | 114713732  | 114713882  | chr1:114.7Mb | 0     | +      | 114713732  | 114713770  |
| chr1       | 114713834  | 114713984  | NRAS         | 0     | -      | 114713948  | 114713984  |
| chr1       | 114713895  | 114714045  | NRAS         | 0     | -      | 114714009  | 114714045  |
| chr1       | 114715969  | 114716119  | NRAS         | 0     | +      | 114715969  | 114716009  |
| chr1       | 114716063  | 114716213  | NRAS         | 0     | -      | 114716183  | 114716213  |

Where Start and End are expected coverage based on read length (in this case, 150 bp), and ThickStart and ThickEnd are primer binding locations. Make a new BED file for just the primer binding locations:

| Chromosome   | Start      | End        | Name     | Score | Strand |
|--------------|------------|------------|----------|-------|--------|
| NC_000001.11 | 114713732  | 114713770  | Primer1  | 20    | +      |
| NC_000001.11 | 114713948  | 114713984  | Primer2  | 20    | -      |
| NC_000001.11 | 114714009  | 114714045  | Primer3  | 20    | -      |
| NC_000001.11 | 114715969  | 114716009  | Primer4  | 20    | +      |
| NC_000001.11 | 114716183  | 114716213  | Primer5  | 20    | -      |

Note chromosome names have been changed to reflect the genome we are using (GCF_GRCh.38).

Using BEDTools, we can get the sequence information with the `getfasta` subcommand.

```
./bedtools getfasta [OPTIONS] -fi <input FASTA> -bed <BED> -fo <output FASTA>
```

We want the output to just contain the sequences of interest, and contain strand information. Therefore, we will use the `-tab` and `-s` options to have the output in a tab-separated format instead of FASTA, and to get sequence from the proper strand:

```
./bedtools getfasta -tab -s -fi <input FASTA> -bed <BED> -fo <output.tsv>
```
Which will give:

| **Genomic Region**                          | **Sequence**                              |
|---------------------------------------------|-------------------------------------------|
| NC_000001.11:114713732-114713770 (+)        | CCTCATTTCCCCATAAAGATTCAGAACACAAAGATCAT    |
| NC_000001.11:114713948-114713984 (-)        | CCCCAGGATTCTTACAGAAAACAAGTGGTTATAGAT      |
| NC_000001.11:114714009-114714045 (-)        | TCCCTgtggtttttaataaaaattgaacttccctcc      |
| NC_000001.11:114715969-114716009 (+)        | CCTTTAATACAGAATATGGGTAAAGATGATCCGACAAGTG  |

This provides the exact primer sequences used during amplification. Note that the QIASeq pro panel only use one primer sequence for target regions at the I7 adapter end (the other is a universal primer matching a ligated sequence), and all reads are formatted 5′ → 3′ based on the I5 adapter. Therefore, the listed primer sequences must be reverse complemented. Masking is also present and must be removed. 

Both steps are performed automatically by the `Cutadapt_primer_trimming.py` script, which accepts this TSV file directly as input.

## Target Variants CSV file
CSV file containing the following headings for each variant of interest. Any variants not listed here will not be examined by ctRCA-Var

| Column Heading                | Description                                                                 |
|------------------------------|-----------------------------------------------------------------------------|
| GENE_NAME                    | Name of the gene where the mutation occurs                                 |
| Mutation CDS                 | Coding DNA sequence (CDS) notation of the mutation (e.g. c.121A>G)                        |
| Mutation AA                  | Amino acid change resulting from the mutation (p.T41A)                              |
| GENOMIC_WT_ALLELE_SEQ        | Reference (wild-type) genomic allele sequence                              |
| GENOMIC_MUT_ALLELE_SEQ       | Mutated genomic allele sequence                                            |
| Mutation Description AA      | Description of the mutation at the amino acid level                        |
| chrom                        | Chromosome on which the mutation is located                                |
| pos                          | Genomic position of the mutation                                           |

The Target mutations info-sheet in the `files` directory provides information on how to how assemble this file if there are many variants