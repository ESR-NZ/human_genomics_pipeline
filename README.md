# human_genomics_pipeline

A Snakemake workflow to process single samples or cohorts of paired-end sequencing data (WGS or WES) using [bwa](http://bio-bwa.sourceforge.net/) and [GATK4](https://gatk.broadinstitute.org/hc/en-us), taking the data from fastq to vcf. Quality control checks are also undertaken. The fastq files can be optionally trimmed and the pipeline can run on [NVIDIA GPU's](https://www.nvidia.com/en-gb/graphics-cards/) where [nvidia clara parabricks software is available](https://www.nvidia.com/en-us/docs/parabricks/quickstart-guide/software-overview/). This workflow is designed to follow the [GATK best practice workflow for germline short variant discovery (SNPs + Indels)](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-) and is designed to be followed by [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline). This pipeline has been developed with human genetic data in mind, however we designed it to be species agnostic. Genetic data from other species can in theory can be analysed by setting a species-specific reference genome and variant databases in the configuration file.

- [human_genomics_pipeline](#human_genomics_pipeline)
  - [Workflow diagram - single samples](#workflow-diagram---single-samples)
  - [Workflow diagram - single samples - GPU accelerated](#workflow-diagram---single-samples---gpu-accelerated)
  - [Workflow diagram - cohort samples](#workflow-diagram---cohort-samples)
  - [Workflow diagram - cohort samples - GPU accelerated](#workflow-diagram---cohort-samples---gpu-accelerated)
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
      - [HPC](#hpc-1)
    - [9. Evaluate the pipeline run](#9-evaluate-the-pipeline-run)
    - [10. Commit and push to your forked version of the github repo](#10-commit-and-push-to-your-forked-version-of-the-github-repo)
    - [11. Repeat step 10 each time you re-run the analysis with different parameters](#11-repeat-step-10-each-time-you-re-run-the-analysis-with-different-parameters)
    - [12. Create a pull request with the upstream repo to merge any useful changes (optional)](#12-create-a-pull-request-with-the-upstream-repo-to-merge-any-useful-changes-optional)

## Workflow diagram - single samples

<img src="./images/rulegraph_single.png" class="center">

## Workflow diagram - single samples - GPU accelerated

<img src="./images/rulegraph_single_gpu.png" class="center">

## Workflow diagram - cohort samples

<img src="./images/rulegraph_cohort.png" class="center">

## Workflow diagram - cohort samples - GPU accelerated

<img src="./images/rulegraph_cohort_gpu.png" class="center">

## Run human_genomics_pipeline

- **Prerequisite hardware:** [NVIDIA GPUs](https://www.nvidia.com/en-gb/graphics-cards/) (for GPU accelerated runs)
- **Prerequisite software:** [NVIDIA CLARA PARABRICKS and dependencies](https://www.nvidia.com/en-us/docs/parabricks/local-installation/) (for GPU accelerated runs), [Git 2.7.4](https://git-scm.com/), [Mamba 0.4.4](https://github.com/TheSnakePit/mamba) with [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html), [gsutil 4.52](https://pypi.org/project/gsutil/), [gunzip 1.6](https://linux.die.net/man/1/gunzip)
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
|     |___proband1_pedigree.ped
|     |___proband2_pedigree.ped
|     |___ ...
|
|___human_genomics_pipeline/

```

Requirements:
  - Input paired end fastq files need to identified with `_1` and `_2` (not `_R1` and `_R2`)
  - Currently, the filenames of the pedigree files need to be labelled with the name of the proband/individual affected with the disease phenotype in the cohort (we will be working towards removing this requirement)
  - Singletons and cohorts need to be run in separate pipeline runs

Assumptions:
  - There is one proband/individual affected with the disease phenotype of interest in a given cohort (one individual with a value of 2 in the 6th column of the pedigree file)

See [here](https://help.github.com/en/github/getting-started-with-github/fork-a-repo#keep-your-fork-synced) for help

### 3. Create a local copy of the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle) (either b37 or hg38)

#### b37

Download from [Google Cloud Bucket](https://console.cloud.google.com/storage/browser/gatk-legacy-bundles/b37?prefix=)

```bash
gsutil cp -r gs://gatk-legacy-bundles/b37/ /where/to/download/
```

#### hg38

Download from [Google Cloud Bucket](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0)

```bash
gsutil cp -r gs://genomics-public-data/resources/broad/hg38/ /where/to/download/
```

### 4. Modify the configuration file

Edit 'config.yaml' found within the config directory.

Specify whether the data is to be analysed on it's own ('Single') or as a part of a cohort of samples ('Cohort'). For example:

```yaml
DATA: "Single"
```

Specify whether the pipeline should be GPU accelerated where possible (either 'Yes' or 'No', this requires [NVIDIA GPUs](https://www.nvidia.com/en-gb/graphics-cards/) and [NVIDIA CLARA PARABRICKS](https://www.nvidia.com/en-us/docs/parabricks/local-installation/))

```yaml
GPU_ACCELERATED: "No"
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

Specify whether the raw fastq reads should be trimmed (either 'Yes' or 'No'). For example:

```yaml
TRIM: "Yes"
```

If trimming the raw fastq reads, set the [trim galore](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md) adapter trimming parameters. Choose one of the common adapters such as Illumina universal, Nextera transposase or Illumina small RNA with `--illumina`, `--nextera` or `--small_rna`. Alternatively, pass adapter sequences to the `-a` and `-a2` flags. If not set, trim galore will try to auto-detect the adapter based on the fastq reads.

```yaml
TRIMMING:
  ADAPTERS: "--illumina"
```

Pass the resources to be used to recalibrate bases with [gatk BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360047217531-BaseRecalibrator) to the `--known-sites` flag if not gpu accelerated and to the `--knownSites` if gpu accelerated. For example:

```yaml
RECALIBRATION:
  RESOURCES: "--known-sites /home/lkemp/publicData/b37/dbsnp_138.b37.vcf
            --known-sites /home/lkemp/publicData/b37/Mills_and_1000G_gold_standard.indels.b37.vcf
            --known-sites /home/lkemp/publicData/b37/1000G_phase1.indels.b37.vcf"
```

### 5. Configure to run on a HPC (optional)

*This will deploy the non-GPU accelerated rules to slurm and deploy the GPU accelerated rules locally (pbrun_fq2bam, pbrun_haplotypecaller_single, pbrun_haplotypecaller_cohort)*

In theory, this cluster configuration should be adaptable to other job scheduler systems, but here I will demonstrate how to deploy this pipeline to [slurm](https://slurm.schedmd.com/).

Configure `account:` and `partition:` in the default section of 'cluster.json' in order to set the parameters for slurm sbatch (see documentation [here](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#cluster-configuration-deprecated) and [here](https://slurm.schedmd.com/)). For example:

```json
{
    "__default__" :
    {
        "account" : "lkemp",
        "nodes" : 1,
        "ntasks" : 4,
        "partition" : "prod"
    },
```

[This](https://hpc-carpentry.github.io/hpc-python/17-cluster/) is a good place to go for a good working example.

### 6. Modify the run scripts

Set the number maximum number of cores to be used with the `-j` flag. If running GPU accelerated, also set the maximum number of GPU's to be used with the `--resources` flag. For example:

Dry run (dryrun.sh):

```bash
snakemake \
-n \
-j 32 \
--resources gpu=2 \
--use-conda \
--conda-frontend mamba \
--configfile ../config/config.yaml
```

Full run (run.sh):

```bash
snakemake \
-j 32 \
--resources gpu=2 \
--use-conda \
--configfile ../config/config.yaml
```

See the [snakemake documentation](https://snakemake.readthedocs.io/en/v4.5.1/executable.html#all-options) for additional run parameters.

#### HPC

If you want to run the pipeline on a HPC, set the `-j` and `--resources` flags in dryrun_hpc.sh and run_hpc.sh run scripts instead. For example:

Dry run (dryrun_hpc.sh):

```bash
snakemake \
-n \
-j 32 \
--resources gpu=2 \
--use-conda \
--conda-frontend mamba \
--configfile ../config/config.yaml \
--cluster-config ../config/cluster.json \
--cluster "sbatch -A {cluster.account} \
-p {cluster.partition} \
--nodes {cluster.nodes} \
--ntasks {cluster.ntasks}"
```

Full run (run_hpc.sh):

```bash
snakemake \
-j 32 \
--resources gpu=2 \
--use-conda \
--conda-frontend mamba \
--configfile ../config/config.yaml \
--cluster-config ../config/cluster.json \
--cluster "sbatch -A {cluster.account} \
-p {cluster.partition} \
--nodes {cluster.nodes} \
--ntasks {cluster.ntasks}"
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

#### HPC

If running on a HPC, run the HPC run scripts instead

```bash
bash dryrun_hpc.sh
bash run_hpc.sh
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
