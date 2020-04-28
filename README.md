# human_genomics_pipeline

A simple Snakemake workflow to process paired-end sequencing data (WGS) using bwa/GATK4.

- [human_genomics_pipeline](#humangenomicspipeline)
  - [workflow diagram](#workflow-diagram)
  - [How to run human_genomics_pipeline](#how-to-run-humangenomicspipeline)
    - [1. Fork the pipeline repo to a personal or lab account](#1-fork-the-pipeline-repo-to-a-personal-or-lab-account)
    - [2. Take the pipeline to the data on your local machine](#2-take-the-pipeline-to-the-data-on-your-local-machine)
    - [3. Create a local copy of the reference human genome and dbSNP database (either GRCh37 or GRCh38)](#3-create-a-local-copy-of-the-reference-human-genome-and-dbsnp-database-either-grch37-or-grch38)
      - [GRCh37](#grch37)
      - [GRCh38](#grch38)
    - [4. Choose and modify an appropriate configuration file](#4-choose-and-modify-an-appropriate-configuration-file)
    - [5. Create and activate a conda environment with python and snakemake and installed](#5-create-and-activate-a-conda-environment-with-python-and-snakemake-and-installed)
    - [6. Run the pipeline](#6-run-the-pipeline)
    - [7. Evaluate the pipeline run](#7-evaluate-the-pipeline-run)
    - [8. Commit and push to your forked version of the repo](#8-commit-and-push-to-your-forked-version-of-the-repo)
    - [9. Create a pull request with the upstream repo to merge any useful changes](#9-create-a-pull-request-with-the-upstream-repo-to-merge-any-useful-changes)
  - [Useful reading](#useful-reading)

## workflow diagram

<img src="rulegraph.png" class="center">

## How to run human_genomics_pipeline

- **Prerequisite software:** [Git](https://git-scm.com/), [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html), [wget](https://www.gnu.org/software/wget/), [gunzip](https://linux.die.net/man/1/gunzip), [bwa](http://bio-bwa.sourceforge.net/)
- **OS:** Validated on Ubuntu 16.04

### 1. Fork the pipeline repo to a personal or lab account

See [here](https://help.github.com/en/github/getting-started-with-github/fork-a-repo#fork-an-example-repository) for help forking a github repository

### 2. Take the pipeline to the data on your local machine

Clone the forked [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) repo into the same directory as your paired end fastq data to be processed. Required folder structure:

```bash

.
|___fastq/
|     |___sample1_R1
|     |___sample1_R2
|     |___sample2_R1
|     |___sample2_R2
|     |___ ...
|
|___human_genomics_pipeline/

```

See [here](https://help.github.com/en/github/getting-started-with-github/fork-a-repo#keep-your-fork-synced) for help cloning a forked github repository.

### 3. Create a local copy of the reference human genome and dbSNP database (either GRCh37 or GRCh38)

#### GRCh37

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

### 4. Choose and modify an appropriate configuration file

- Use 'config_GRCh37.yaml' to run the pipeline against the GRCh37 reference genome
- Use 'config_GRCh38.yaml' to run the pipeline against the GRCh38 reference genome

Set the the working directories in the config file to the reference human genome file, dbSNP database file and a temporary directory. For example:

```yaml
GENOME: "/home/lkemp/publicData/referenceGenome/Homo_sapiens_assembly38.fasta.gz"
dbSNP: "/home/lkemp/publicData/dbSNP/All_20180418.vcf.gz"
TEMPDIR: "/home/lkemp/tmp/"
```

### 5. Create and activate a conda environment with python and snakemake and installed

```bash
conda create -n pipeline_env python=3.7
conda activate pipeline_env
conda install -c bioconda snakemake=5.14.0
```

### 6. Run the pipeline

Specify the config file to be used with the `--configfile` flag and modify the number of cores to be used with the `-j` flag. First carry out a dry run. If there are no issues, start a full run without the `-n` flag.

Dry run:

```bash
snakemake -n -j 24 --use-conda --configfile config.yaml
```

Full run:

```bash
snakemake -j 24 --use-conda --configfile config.yaml
```

See the [snakemake documentation](https://snakemake.readthedocs.io/en/v4.5.1/executable.html) for additional run parameters.

### 7. Evaluate the pipeline run

```bash
snakemake --report report.html --configfile config.yaml
```

### 8. Commit and push to your forked version of the repo

To maintain reproducibility, commit and push:

- All modified configuration file/s
- Output from the run such as reports and plots (optional)

### 9. Create a pull request with the [upstream repo](https://github.com/ESR-NZ/human_genomics_pipeline) to merge any useful changes

Contributions and feedback are more than welcome! :blush:

See [here](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request) for help creating a pull request.

## Useful reading

Van der Auwera et al., (2013). *Current Protocols in Bioinformatics*. [From FastQ data to high confidence variant calls: the Genome Analysis Toolkit best practices pipeline](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4243306/); 11(1110): 11.10.1–11.10.33.
