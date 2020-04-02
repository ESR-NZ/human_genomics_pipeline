# human_genomics_pipeline

A simple Snakemake workflow to process paired-end sequencing data (WGS) using bwa/GATK4.

## Table of contents

- [human_genomics_pipeline](#humangenomicspipeline)
  - [Table of contents](#table-of-contents)
  - [workflow diagram](#workflow-diagram)
  - [Set up and run vcf_annotation_pipeline](#set-up-and-run-vcfannotationpipeline)
    - [Download pipeline, reference genome and dbSNP](#download-pipeline-reference-genome-and-dbsnp)
      - [GRCh37](#grch37)
      - [GRCh38](#grch38)
    - [Set up the working environment](#set-up-the-working-environment)
    - [Run the pipeline](#run-the-pipeline)
    - [Resource allocation](#resource-allocation)
  - [Useful links/papers](#useful-linkspapers)

## workflow diagram

<img src="rulegraph.png" class="center">

## Set up and run vcf_annotation_pipeline

- **Prerequisite software:** [Git](https://git-scm.com/), [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html), [gunzip](https://linux.die.net/man/1/gunzip), [bwa](http://bio-bwa.sourceforge.net/)
- **Prerequisite data:** None
- **OS:** Validated on Ubuntu 16.04

### Download pipeline, reference genome and dbSNP

Clone the [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) repository

```bash
git clone https://github.com/ESR-NZ/human_genomics_pipeline.git
```

#### GRCh37

Create a GRCh37 folder within a publicData folder to download the reference genome and dbSNP into

```bash
mkdir publicData
cd publicData
mkdir GRCh37
cd GRCh37
```

Download the reference human genome (GRCh37) and it's associated fasta sequence dictionary file (.dict) and fasta index file (.fai) files from the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)

```bash
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.fasta.gz
gunzip ucsc.hg19.fasta.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.dict.gz
gunzip ucsc.hg19.dict.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.fasta.fai.gz
gunzip ucsc.hg19.fasta.fai.gz
```

Create index files for the genome sequence (.amb, .ann, .bwt, .pac, .sa)

```bash
bwa index -a bwtsw ucsc.hg19.fasta
```

Download dbSNP (build 151) from [NCBI](https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/)

```bash
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/GATK/All_20180423.vcf.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/GATK/All_20180423.vcf.gz.tbi
```

#### GRCh38

Create a GRCh38 folder within a publicData folder to download the reference genome and dbSNP into

```bash
mkdir publicData
cd publicData
mkdir GRCh38
cd GRCh38
```

Download the reference human genome (GRCh38) and it's associated fasta sequence dictionary file (.dict) and fasta index file (.fai) files from the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)

```bash
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.gz
gunzip Homo_sapiens_assembly38.fasta.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.dict
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.fai
```

Create index files for the genome sequence (.amb, .ann, .bwt, .pac, .sa)

```bash
bwa index -a bwtsw Homo_sapiens_assembly38.fasta
```

Download dbSNP (build 151) from [NCBI](https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/)

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/All_20180418.vcf.gz
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/All_20180418.vcf.gz.tbi
```

### Set up the working environment

Set the the appropriate variables in 'config.yaml'. Choose to run the pipeline against either GRCh37 or GRCh38 by setting the BUILD variable

```yaml
BUILD:
  "GRCh38"
```

Also set the directory to your temporary files folder

```yaml
TEMPDIR:
  "/store/lkemp/exome_project/tmp/"
```

And the file directories to the publicData folder containing the data we downloaded above

```yaml
PUBLICDIR:
  "/store/lkemp/publicData/"
```

Lastly, set the file directory to the sequence data (fastq files)

```yaml
SAMPLEDIR:
  "/store/lkemp/exome_project/data/fastq/"
```

Create and activate a conda environment with python and snakemake and installed

```bash
conda create -n pipeline_env python=3.7
conda activate pipeline_env
conda install -c bioconda snakemake=5.10.0
```

### Run the pipeline

First start a dry run. If there are no issues, start a full run without the -n flag

```bash
snakemake -n -r -j 24 -p --use-conda
snakemake -r -j 24 -p --use-conda
```

### Resource allocation

If necessary, the maximum number of CPU cores allocated by changing the -j flag in the snakemake program. For example to scale to run on a laptop/desktop...

```bash
snakemake -r -j 4 -p --use-conda --use-singularity
```

## Useful links/papers

Van der Auwera et al., (2013). *Current Protocols in Bioinformatics*. [From FastQ data to high confidence variant calls: the Genome Analysis Toolkit best practices pipeline](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4243306/); 11(1110): 11.10.1–11.10.33.
