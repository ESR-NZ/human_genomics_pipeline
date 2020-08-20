# human_genomics_pipeline

A simple Snakemake workflow to process paired-end sequencing data (WGS or WES) using [bwa](http://bio-bwa.sourceforge.net/) and [GATK4](https://gatk.broadinstitute.org/hc/en-us).

- [human_genomics_pipeline](#human_genomics_pipeline)
  - [Workflow diagram - single samples](#workflow-diagram---single-samples)
  - [Workflow diagram - cohorts](#workflow-diagram---cohorts)
  - [Run human_genomics_pipeline](#run-human_genomics_pipeline)
    - [1. Fork the pipeline repo to a personal or lab account](#1-fork-the-pipeline-repo-to-a-personal-or-lab-account)
    - [2. Take the pipeline to the data on your local machine](#2-take-the-pipeline-to-the-data-on-your-local-machine)
    - [3. Create a local copy of the GATK resource bundle (either b37 or hg38)](#3-create-a-local-copy-of-the-gatk-resource-bundle-either-b37-or-hg38)
      - [b37](#b37)
      - [hg38](#hg38)
    - [4. Modify the configuration file](#4-modify-the-configuration-file)
    - [5. Configure to run on a HPC (optional)](#5-configure-to-run-on-a-hpc-optional)
    - [6. Modify the run scripts](#6-modify-the-run-scripts)
      - [HPC](#hpc)
    - [7. Create and activate a conda environment with python and snakemake installed](#7-create-and-activate-a-conda-environment-with-python-and-snakemake-installed)
    - [8. Run the pipeline](#8-run-the-pipeline)
    - [9. Evaluate the pipeline run](#9-evaluate-the-pipeline-run)
    - [10. Commit and push to your forked version of the github repo](#10-commit-and-push-to-your-forked-version-of-the-github-repo)
    - [11. Repeat step 10 each time you re-run the analysis with different parameters](#11-repeat-step-10-each-time-you-re-run-the-analysis-with-different-parameters)
    - [12. Create a pull request with the upstream repo to merge any useful changes (optional)](#12-create-a-pull-request-with-the-upstream-repo-to-merge-any-useful-changes-optional)

## Workflow diagram - single samples

<img src="./images/rulegraph_single.png" class="center">

## Workflow diagram - cohorts

<img src="./images/rulegraph_cohort.png" class="center">

## Run human_genomics_pipeline

- **Prerequisite software:** [Git 2.7.4](https://git-scm.com/), [Mamba 0.4.4](https://github.com/TheSnakePit/mamba) with [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html), [gsutil 4.52](https://pypi.org/project/gsutil/), [gunzip 1.6](https://linux.die.net/man/1/gunzip)
- **OS:** Validated on Ubuntu 16.04

### 1. Fork the pipeline repo to a personal or lab account

See [here](https://help.github.com/en/github/getting-started-with-github/fork-a-repo#fork-an-example-repository) for help

### 2. Take the pipeline to the data on your local machine

Clone the forked [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) repo into the same directory as your paired end fastq data to be processed. Required folder structure and file naming convention:

```bash

.
|___fastq/
|     |___sample1_1.fastq.gz
|     |___sample1_2.fastq.gz
|     |___sample2_1.fastq.gz
|     |___sample2_2.fastq.gz
|     |___ ...
|
|___human_genomics_pipeline/

```

If you are analysing cohort's of samples, you will need an additional directory with a [pedigree file](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format) for each cohort/family:

```bash

.
|___fastq/
|     |___sample1_1.fastq.gz
|     |___sample1_2.fastq.gz
|     |___sample2_1.fastq.gz
|     |___sample2_2.fastq.gz
|     |___ ...
|
|_____pedigrees/
|     |___family1_pedigree.ped
|     |___family2_pedigree.ped
|     |___ ...
|
|___human_genomics_pipeline/

```

Requirements:
  - Input paired end fastq files need to identified with `_1` and `_2` (not `_R1` and `_R2`)
  - The filenames of the pedigree files need to be labelled with the family name
  - Singletons and cohorts need to be run in seperate pipeline runs

See [here](https://help.github.com/en/github/getting-started-with-github/fork-a-repo#keep-your-fork-synced) for help

### 3. Create a local copy of the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle) (either b37 or hg38)

#### b37

Download from [Google Cloud Bucket](https://console.cloud.google.com/storage/browser/gatk-legacy-bundles/b37?prefix=)

```bash
gsutil cp -r gs://gatk-legacy-bundles/b37 /where/to/download/
```

Unzip all zipped files

```bash
gunzip -f /location/you/downloaded/bundle/*.gz
```

#### hg38

Download from [Google Cloud Bucket](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0)

```bash
gsutil cp -r gs://genomics-public-data/resources/broad/hg38/v0/ /where/to/download/
```

Unzip all zipped files

```bash
gunzip -f /location/you/downloaded/bundle/*.gz
```

### 4. Modify the configuration file

Edit 'config.yaml' found within the config directory.

Specify whether the data is to be analysed on it's own ('Single') or as a part of a cohort of samples ('Cohort'). For example:

```yaml
DATA: "Single"
```

Set the the working directories to the reference human genome file (b37 or hg38) and it's associated dictionary file (.dict). For example:

```yaml
REFGENOME: "/home/lkemp/publicData/b37/human_g1k_v37_decoy.fasta"
DICT: "/home/lkemp/publicData/b37/human_g1k_v37_decoy.dict"
```

Set the the working directory to your dbSNP database file (b37 or hg38). For example:

```yaml
dbSNP: "/home/lkemp/publicData/b37/dbsnp_138.b37.vcf"
```

Set the the working directory to a temporary file directory. For example:

```yaml
TEMPDIR: "/home/lkemp/tmp/"
```

If analysing WES data, pass a design file (.bed) indicating the genomic regions that were sequenced to the `-L` flag (see [here](https://leahkemp.github.io/documentation/human_genomic_pipelines/design_files.html) for more information on accessing design files). Also set the level of padding by passing the amount of padding in base pairs to the `-ip` flag. For example:

*If NOT analysing WES data, leave these fields blank*

```yaml
WES:
  # File path to the exome capture regions over which to operate (prefix with the '-L' flag)
  INTERVALS: "-L /home/lkemp/publicData/sure_select_human_all_exon_V7/S31285117_Padded.bed"
  # Padding (in bp) to add to each region (prefix with the '-ip' flag)
  PADDING: "-ip 100"
```

Pass the resources to be used to recalibrate bases with [gatk BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360047217531-BaseRecalibrator) to the `--known-sites` flag. For example:

```yaml
RECALIBRATION:
  RESOURCES: "--known-sites /home/lkemp/publicData/b37/dbsnp_138.b37.vcf
            --known-sites /home/lkemp/publicData/b37/Mills_and_1000G_gold_standard.indels.b37.vcf
            --known-sites /home/lkemp/publicData/b37/1000G_phase1.indels.b37.vcf"
```

### 5. Configure to run on a HPC (optional)

In theory, this cluster configuration should be adaptable to other job scheduler systems, but here I will demonstrate how to deploy this pipeline to [slurm](https://slurm.schedmd.com/).

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

[This](https://hpc-carpentry.github.io/hpc-python/17-cluster/) is a good place to go for a good working example.

*These variables will need to be passed to snakemake in the snakemake run script (see example in the HPC section of step 6).*

### 6. Modify the run scripts

Set the number of cores to be used with the `-j` flag. For example:

Dry run (dryrun.sh):

```bash
snakemake \
-n \
-j 32 \
--use-conda \
--configfile ../config/config.yaml
```

Full run (run.sh):

```bash
snakemake \
-j 32 \
--use-conda \
--configfile ../config/config.yaml
```

Report (report.sh)

```bash
snakemake \
--report ../results/report.html \
--configfile ../config/config.yaml \
--report-stylesheet ../config/ESR_stylesheet.css
```

See the [snakemake documentation](https://snakemake.readthedocs.io/en/v4.5.1/executable.html#all-options) for additional run parameters.

#### HPC

If you want to run the pipeline on a HPC, pass the cluster variables set in 'cluster.json' to the dry run and full run scripts. For example:

Dry run (dryrun.sh):

```bash
snakemake \
-n \
-j 32 \
--use-conda \
--configfile ../config/config.yaml \
--cluster-config ../config/cluster.json \
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
--configfile ../config/config.yaml \
--cluster-config ../config/cluster.json \
--cluster "sbatch -A {cluster.account} \
-p {cluster.partition} \
--nodes {cluster.nodes} \
--ntasks {cluster.ntasks} \
--nodelist {cluster.nodelist}"
```

### 7. Create and activate a conda environment with python and snakemake installed

```bash
cd ./human_genomics_pipeline/workflow/
mamba env create -f pipeline_run_env.yml
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

### 10. Commit and push to your forked version of the github repo

To maintain reproducibility, commit and push:

- All configuration files
- All run scripts
- The final report

### 11. Repeat step 10 each time you re-run the analysis with different parameters

### 12. Create a pull request with the [upstream repo](https://github.com/ESR-NZ/human_genomics_pipeline) to merge any useful changes (optional)

Contributions and feedback are more than welcome! :blush:

See [here](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request) for help
