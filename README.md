# human_genomics_pipeline

A simple Snakemake workflow to process paired-end sequencing data (WGS or WES) using [bwa](http://bio-bwa.sourceforge.net/), [sambamba](https://lomereiter.github.io/sambamba/) and [GATK4](https://gatk.broadinstitute.org/hc/en-us).

- [human_genomics_pipeline](#human_genomics_pipeline)
  - [workflow diagram](#workflow-diagram)
  - [How to run human_genomics_pipeline](#how-to-run-human_genomics_pipeline)
    - [1. Fork the pipeline repo to a personal or lab account](#1-fork-the-pipeline-repo-to-a-personal-or-lab-account)
    - [2. Take the pipeline to the data on your local machine](#2-take-the-pipeline-to-the-data-on-your-local-machine)
    - [3. Create a local copy of the reference human genome and variation databases (either GRCh37 or GRCh38)](#3-create-a-local-copy-of-the-reference-human-genome-and-variation-databases-either-grch37-or-grch38)
      - [GRCh37](#grch37)
      - [GRCh38](#grch38)
    - [4. Modify the configuration file](#4-modify-the-configuration-file)
    - [5. Configure to run on a HPC (optional)](#5-configure-to-run-on-a-hpc-optional)
    - [6. Modify the run scripts](#6-modify-the-run-scripts)
      - [HPC](#hpc)
    - [7. Create and activate a conda environment with python and snakemake installed](#7-create-and-activate-a-conda-environment-with-python-and-snakemake-installed)
    - [8. Run the pipeline](#8-run-the-pipeline)
    - [9. Evaluate the pipeline run](#9-evaluate-the-pipeline-run)
    - [10. Commit and push to your forked version of the repo](#10-commit-and-push-to-your-forked-version-of-the-repo)
    - [11. Repeat step 10 each time you re-run the analysis with different parameters](#11-repeat-step-10-each-time-you-re-run-the-analysis-with-different-parameters)
    - [12. Create a pull request with the upstream repo to merge any useful changes (optional)](#12-create-a-pull-request-with-the-upstream-repo-to-merge-any-useful-changes-optional)
  - [Useful reading](#useful-reading)

## workflow diagram

<img src="rulegraph.png" class="center">

## How to run human_genomics_pipeline

- **Prerequisite software:** [Git](https://git-scm.com/), [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html), [wget](https://www.gnu.org/software/wget/), [gunzip](https://linux.die.net/man/1/gunzip), [bwa](http://bio-bwa.sourceforge.net/)
- **OS:** Validated on Ubuntu 16.04

---

*This pipeline is in ongoing development. It is currently operational for single samples of whole exome sequencing data (not cohort data yet). It should theoretically be operational for whole genome sequencing data, but this has not yet been validated. It can be run against either the GRCh37 or GRCh38 reference genome. Benchmarking is currently underway.*

---

### 1. Fork the pipeline repo to a personal or lab account

See [here](https://help.github.com/en/github/getting-started-with-github/fork-a-repo#fork-an-example-repository) for help

### 2. Take the pipeline to the data on your local machine

Clone the forked [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) repo into the same directory as your paired end fastq data to be processed. Required folder structure and file naming convention:

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

See [here](https://help.github.com/en/github/getting-started-with-github/fork-a-repo#keep-your-fork-synced) for help

### 3. Create a local copy of the reference human genome and variation databases (either GRCh37 or GRCh38)

#### GRCh37

Download the reference human genome (GRCh37) and it's associated fasta sequence dictionary file (.dict) and fasta index file (.fai) files from the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)

```bash
wget ftp://gsapubftp-anonymous:password@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.fasta.gz
gunzip ucsc.hg19.fasta.gz
wget ftp://gsapubftp-anonymous:password@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.dict.gz
gunzip ucsc.hg19.dict.gz
wget ftp://gsapubftp-anonymous:password@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.fasta.fai.gz
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

Download other variation databases from the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)

```bash
wget ftp://gsapubftp-anonymous:password@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
wget ftp://gsapubftp-anonymous:password@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx.gz
wget ftp://gsapubftp-anonymous:password@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.gz
wget ftp://gsapubftp-anonymous:password@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.idx.gz
```

#### GRCh38

Download the reference human genome (GRCh38) and it's associated fasta sequence dictionary file (.dict) and fasta index file (.fai) files from the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)

```bash
wget ftp://gsapubftp-anonymous:password@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.gz
gunzip Homo_sapiens_assembly38.fasta.gz
wget ftp://gsapubftp-anonymous:password@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.dict
wget ftp://gsapubftp-anonymous:password@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.fai
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

Download other variation databases from the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)

```bash
wget ftp://gsapubftp-anonymous:password@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget ftp://gsapubftp-anonymous:password@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
```

### 4. Modify the configuration file

Either edit the 'config_template.yaml' file from scratch, or edit a config file 'config_examples'.

Specify whether you are running your analysis against the GRCh37 or GRCh38 build of the reference genome and whether the data is to be analysed on it's own ('Single') or as a part of a cohort ('Cohort'). For example:

```yaml
# Type of input data (either 'Single' or 'Cohort')
DATA: "Single"
```

Set the the working directories in the config file to the reference human genome file (GRCh37 or GRCh38), dbSNP database file (GRCh37 or GRCh38) and a temporary directory. For example:

```yaml
# File directories to reference genome and dbSNP database
REFGENOME: "/home/lkemp/publicData/referenceGenome/Homo_sapiens_assembly38.fasta"
dbSNP: "/home/lkemp/publicData/dbSNP/All_20180418.vcf.gz"

# Temporary file directory
TEMPDIR: "/home/lkemp/tmp/"
```

If analysing WES data, pass a design file (.bed) indicating the genomic regions that were sequenced (see [here](https://leahkemp.github.io/documentation/human_genomic_pipelines/design_files.html) for more information on accessing design files). Also set the level of padding. For example:

*If NOT analysing WES data, leave these fields blank*

```yaml
# WES settings (leave blank if analysing other data such as WGS)
WES:
  # Genomic intervals over which to operate
  INTERVALS: "-L /home/lkemp/publicData/sure_select_human_all_exon_V7/S31285117_AllTracks.bed"
  # Padding (in bp) to add to each interval
  PADDING: "-ip 100"
```

Pass the resources to be used to recalibrate bases with gatk BaseRecalibrator to the `--known-sites` flag. For example:

```yaml
RECALIBRATION:
  RESOURCES: "--known-sites /home/lkemp/publicData/dbSNP/All_20180418.vcf.gz
              --known-sites /home/lkemp/publicData/mills/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
```

Save your modified config file with a descriptive name

### 5. Configure to run on a HPC (optional)

In theory, this cluster configuration should be adaptable to other job scheduler systems, but here I will demonstrate how to integrate this pipeline with slurm.

Configure `account:`, `partition:` and `nodelist:` in default section of 'cluster.json' in order to set the parameters for slurm sbatch (see documentation [here](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#cluster-configuration-deprecated) and [here](https://slurm.schedmd.com/)). For example:

```json
{
    "__default__" :
    {
        "account" : "lkemp",
        "nodes" : 1,
        "ntasks" : 4,
        "partition" : "prod",
        "nodelist" : "kscprod-bio4"
    }
}
```

[This](https://hpc-carpentry.github.io/hpc-python/17-cluster/) is a good place to go for a good working example if you get stuck.

*These variables will need to be passed to snakemake in the snakemake run script (see example in step 6).*

### 6. Modify the run scripts

Specify your config file to be used with the `--configfile` flag and modify the number of cores to be used with the `-j` flag. For example:

Dry run (dryrun.sh):

```bash
snakemake \
-n -j 32 \
--use-conda \
--configfile config_37_single_WES.yaml
```

Full run (run.sh):

```bash
snakemake \
-j 32 \
--use-conda \
--configfile config_37_single_WES.yaml
```

Report (report.sh)

```bash
snakemake \
--report ESR_report.html \
--configfile config_37_single_WES.yaml \
--report-stylesheet ESR_stylesheet.css
```

See the [snakemake documentation](https://snakemake.readthedocs.io/en/v4.5.1/executable.html) for additional run parameters.

#### HPC

If you want to run the pipeline on a HPC, pass the cluster variables set in 'cluster.json' to the dry run and full run scripts. For example:

Dry run (dryrun.sh):

```bash
snakemake \
-n -j 32 \
--use-conda \
--configfile config_37_single_WES.yaml \
--cluster-config cluster.json \
--cluster "sbatch -A {cluster.account} \
-p {cluster.partition} \
--nodes {cluster.nodes} \
--ntasks {cluster.ntasks} \
--nodelist {cluster.nodelist}"
```

Full run (run.sh):

```bash
snakemake \
-j 32 \
--use-conda \
--configfile config_37_single_WES.yaml \
--cluster-config cluster.json \
--cluster "sbatch -A {cluster.account} \
-p {cluster.partition} \
--nodes {cluster.nodes} \
--ntasks {cluster.ntasks} \
--nodelist {cluster.nodelist}"
```

### 7. Create and activate a conda environment with python and snakemake installed

```bash
conda env create -f pipeline_run_env.yml
conda activate pipeline_run_env
```

### 8. Run the pipeline

First carry out a dry run

```bash
bash dryrun.sh
```

If there are no issues, start a full run

```bash
bash run.sh
```

### 9. Evaluate the pipeline run

Generate an interactive html report

```bash
bash report.sh
```

### 10. Commit and push to your forked version of the repo

To maintain reproducibility, commit and push:

- All configuration files
- All run scripts
- The final report

### 11. Repeat step 10 each time you re-run the analysis with different parameters

### 12. Create a pull request with the [upstream repo](https://github.com/ESR-NZ/human_genomics_pipeline) to merge any useful changes (optional)

Contributions and feedback are more than welcome! :blush:

See [here](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request) for help

## Useful reading

Van der Auwera et al., (2013). *Current Protocols in Bioinformatics*. [From FastQ data to high confidence variant calls: the Genome Analysis Toolkit best practices pipeline](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4243306/); 11(1110): 11.10.1–11.10.33.
