# Run human_genomics_pipeline on a single machine like a laptop or single server/computer

## Table of contents

- [Run human_genomics_pipeline on a single machine like a laptop or single server/computer](#run-human_genomics_pipeline-on-a-single-machine-like-a-laptop-or-single-servercomputer)
  - [Table of contents](#table-of-contents)
  - [1. Fork the pipeline repo to a personal or lab account](#1-fork-the-pipeline-repo-to-a-personal-or-lab-account)
  - [2. Take the pipeline to the data on your local machine](#2-take-the-pipeline-to-the-data-on-your-local-machine)
  - [3. Setup files and directories](#3-setup-files-and-directories)
    - [Test data](#test-data)
  - [4. Get prerequisite software/hardware](#4-get-prerequisite-softwarehardware)
  - [5. Create and activate a conda environment with software for downloading databases](#5-create-and-activate-a-conda-environment-with-software-for-downloading-databases)
  - [6. Create a local copy of the GATK resource bundle (either b37 or hg38)](#6-create-a-local-copy-of-the-gatk-resource-bundle-either-b37-or-hg38)
    - [b37](#b37)
    - [hg38](#hg38)
  - [7. Modify the configuration file](#7-modify-the-configuration-file)
    - [Overall workflow](#overall-workflow)
    - [Pipeline resources](#pipeline-resources)
      - [Trimming](#trimming)
    - [Base recalibration](#base-recalibration)
  - [8. Modify the run scripts](#8-modify-the-run-scripts)
  - [9. Create and activate a conda environment with software for running the pipeline](#9-create-and-activate-a-conda-environment-with-software-for-running-the-pipeline)
  - [10. Run the pipeline](#10-run-the-pipeline)
  - [11. Evaluate the pipeline run](#11-evaluate-the-pipeline-run)
  - [12. Commit and push to your forked version of the github repo](#12-commit-and-push-to-your-forked-version-of-the-github-repo)
  - [13. Repeat step 12 each time you re-run the analysis with different parameters](#13-repeat-step-12-each-time-you-re-run-the-analysis-with-different-parameters)
  - [14. Raise issues, create feature requests or create a pull request with the upstream repo to merge any useful changes to the pipeline (optional)](#14-raise-issues-create-feature-requests-or-create-a-pull-request-with-the-upstream-repo-to-merge-any-useful-changes-to-the-pipeline-optional)

## 1. Fork the pipeline repo to a personal or lab account

See [here](https://help.github.com/en/github/getting-started-with-github/fork-a-repo#fork-an-example-repository) for help forking a repository

## 2. Take the pipeline to the data on your local machine

Clone the forked [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) repo into the same directory as your paired end fastq data to be processed.

See [here](https://help.github.com/en/github/getting-started-with-github/fork-a-repo#keep-your-fork-synced) for help cloning a repository

## 3. Setup files and directories

Required folder structure and file naming convention:

```bash

.
|___fastq/
|     |___sample1_1.fastq.gz
|     |___sample1_2.fastq.gz
|     |___sample2_1.fastq.gz
|     |___sample2_2.fastq.gz
|     |___sample3_1.fastq.gz
|     |___sample3_2.fastq.gz
|     |___sample4_1.fastq.gz
|     |___sample4_2.fastq.gz
|     |___sample5_1.fastq.gz
|     |___sample5_2.fastq.gz
|     |___sample6_1.fastq.gz
|     |___sample6_2.fastq.gz
|     |___ ...
|
|___human_genomics_pipeline/

```

If you're analysing cohort's of samples, you will need an additional directory with a [pedigree file](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format) for each cohort/family using the following folder structure and file naming convention:

```bash

.
|___fastq/
|     |___sample1_1.fastq.gz
|     |___sample1_2.fastq.gz
|     |___sample2_1.fastq.gz
|     |___sample2_2.fastq.gz
|     |___sample3_1.fastq.gz
|     |___sample3_2.fastq.gz
|     |___sample4_1.fastq.gz
|     |___sample4_2.fastq.gz
|     |___sample5_1.fastq.gz
|     |___sample5_2.fastq.gz
|     |___sample6_1.fastq.gz
|     |___sample6_2.fastq.gz
|     |___ ...
|
|___pedigrees/
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
- It is assumed that there is one proband/individual affected with the disease phenotype of interest in a given cohort (one individual with a value of 2 in the 6th column of a given pedigree file)

### Test data

The provided [test dataset](./test) can be used. Setup the test dataset before running the pipeline on this data - choose to setup to run either a single sample analysis or a cohort analysis with the `-a` flag. For example:

```bash
cd ./human_genomics_pipeline
bash ./test/setup_test.sh -a cohort
```

## 4. Get prerequisite software/hardware

For GPU accelerated runs, you'll need [NVIDIA GPUs](https://www.nvidia.com/en-gb/graphics-cards/) (tested with NVIDIA V100) and [NVIDIA CLARA PARABRICKS and dependencies](https://www.nvidia.com/en-us/docs/parabricks/local-installation/) (tested with parabricks version 3.6.1-1). Talk to your system administrator to see if the HPC has this hardware and software available.

Other software required to get setup and run the pipeline:

- [Git](https://git-scm.com/) (tested with version 1.8.3.1)
- [Conda](https://docs.conda.io/projects/conda/en/latest/index.html) (tested with version 4.11.0)
- [Mamba](https://github.com/TheSnakePit/mamba) (tested with version 0.19.1) (note. [mamba can be installed via conda with a single command](https://mamba.readthedocs.io/en/latest/installation.html#existing-conda-install))

This software is commonly pre-installed on HPC's, likely available as modules that can be loaded. Talk to your system administrator if you need help with this.

## 5. Create and activate a conda environment with software for downloading databases

This installs [gsutil](https://cloud.google.com/storage/docs/gsutil) and it's dependencies

```bash
cd ./workflow/
mamba env create -f ./envs/hgp_download_db_env.yaml
conda activate hgp_download_db_env
```

## 6. Create a local copy of the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle) (either b37 or hg38)

### b37

Download from [Google Cloud Bucket](https://console.cloud.google.com/storage/browser/gatk-legacy-bundles/b37?prefix=)

```bash
gsutil -m cp -r gs://gatk-legacy-bundles/b37 /where/to/download/
```

### hg38

Download from [Google Cloud Bucket](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0)

```bash
gsutil -m cp -r gs://genomics-public-data/resources/broad/hg38 /where/to/download/
```

## 7. Modify the configuration file

Edit 'config.yaml' found within the config directory

### Overall workflow

Specify whether the data is to be analysed on it's own ('Single') or as a part of a cohort of samples ('Cohort'). For example:

```yaml
DATA: "Single"
```

Specify whether the pipeline should be GPU accelerated where possible (either 'Yes' or 'No', this requires [NVIDIA GPUs](https://www.nvidia.com/en-gb/graphics-cards/) and [NVIDIA CLARA PARABRICKS](https://www.nvidia.com/en-us/docs/parabricks/local-installation/))

```yaml
GPU_ACCELERATED: "No"
```

Set the the working directories to the reference human genome file (b37 or hg38). For example:

```yaml
REFGENOME: "/scratch/publicData/b37/human_g1k_v37_decoy.fasta"
```

Set the the working directory to your dbSNP database file (b37 or hg38). For example:

```yaml
dbSNP: "/scratch/publicData/b37/dbsnp_138.b37.vcf"
```

Set the the working directory to a temporary file directory. Make sure this is a location with a fair amount of memory space for large intermediate analysis files. For example:

```yaml
TEMPDIR: "/scratch/tmp/"
```

If analysing WES data, pass a design file (.bed) indicating the genomic regions that were sequenced (see [here](https://leahkemp.github.io/documentation/human_genomic_pipelines/design_files.html) for more information on accessing design files). Also set the level of padding by passing the amount of padding in base pairs. For example:

*If NOT analysing WES data, leave these fields blank*

```yaml
WES:
  # File path to the exome capture regions over which to operate
  INTERVALS: "/scratch/publicData/sure_select_human_all_exon_V7/S31285117_Padded.bed"
  # Padding (in bp) to add to each region
  PADDING: "100"
```

### Pipeline resources

These settings allow you to configure the resources per rule/sample

Set the number of threads to use per sample/rule for multithreaded rules (`rule trim_galore_pe` and `rule bwa_mem`). Multithreading will significantly speed up these rules, however the improvements in speed will diminish beyond 8 threads. If desired, a different number of threads can be set for these multithreaded rules by utilising the `--set-threads` flag in the runscript (see [step 8](#8-modify-the-run-scripts)).

```yaml
THREADS: 8
```

Set the maximum memory usage per rule/sample (eg. '40g' for 40 gigabytes, this should suffice for exomes)

```yaml
MAXMEMORY: "40g"
```

Set the maximum number of GPU's to be used per rule/sample for gpu-accelerated runs (eg `1` for 1 GPU)

```yaml
GPU: 1
```

It is a good idea to consider the number of samples that you are processing. For example, if you set `THREADS: "8"` and set the maximum number of jobs to be run in the run script to `-j/--jobs 32` (see [step 8](#8-modify-the-run-scripts)), a maximum of 3 samples will be able to run at one time for these rules (if they are deployed at the same time), but each sample will complete faster. In contrast, if you set `THREADS: "1"` and `-j/--jobs 32`, a maximum of 32 samples could be run at one time, but each sample will take longer to complete. This also needs to be considered when setting `MAXMEMORY` + `--resources mem_mb` and `GPU` + `--resources gpu`.

#### Trimming

Specify whether the raw fastq reads should be trimmed (either 'Yes' or 'No'). For example:

```yaml
TRIM: "Yes"
```

If trimming the raw fastq reads, set the [trim galore](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md) adapter trimming parameters. Choose one of the common adapters such as Illumina universal, Nextera transposase or Illumina small RNA with `--illumina`, `--nextera` or `--small_rna`. Alternatively, pass adapter sequences to the `-a` and `-a2` flags. If not set, trim galore will try to auto-detect the adapter based on the fastq reads.

```yaml
TRIMMING:
  ADAPTERS: "--illumina"
```

### Base recalibration

Pass the resources to be used to recalibrate bases with [gatk BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360036726891-BaseRecalibrator) (this passes the resource files to the [`--known-sites`](https://gatk.broadinstitute.org/hc/en-us/articles/360036726891-BaseRecalibrator#--known-sites) flag), these known polymorphic sites will be used to exclude regions around known polymorphisms from base recalibration. Note. you can include as many or as few resources as you like, but you'll need at least one recalibration resource file. For example:

```yaml
RECALIBRATION:
  RESOURCES:
    - /scratch/publicData/b37/dbsnp_138.b37.vcf
    - /scratch/publicData/b37/Mills_and_1000G_gold_standard.indels.b37.vcf
    - /scratch/publicData/b37/1000G_phase1.indels.b37.vcf
```

## 8. Modify the run scripts

Set the singularity bind location to a directory that contains your pipeline working directory with the `--singularity-args '-B'` flag. Set the number maximum number of job to be deployed with the `--jobs` flag and the maximum amount of memory to be used (in megabytes) with the `resources mem_mb=` flag. If running GPU accelerated, also set the maximum number of GPU's to be used with the `--resources gpu=` flag. For example:

Dry run (dryrun.sh):

```bash
#!/bin/bash -x

snakemake \
--dryrun \
--jobs 32 \
--resources mem_mb=150000 \
--resources gpu=2 \
--use-conda \
--conda-frontend mamba \
--latency-wait 120 \
--use-singularity \
--singularity-args '-B /bind/location/' \
--configfile ../config/config.yaml
```

Full run (run.sh):

```bash
#!/bin/bash -x

snakemake \
--jobs 32 \
--resources mem_mb=150000 \
--resources gpu=2 \
--use-conda \
--conda-frontend mamba \
--latency-wait 120 \
--use-singularity \
--singularity-args '-B /bind/location/' \
--configfile ../config/config.yaml
```

See the [snakemake documentation](https://snakemake.readthedocs.io/en/v4.5.1/executable.html#all-options) for additional run parameters.

## 9. Create and activate a conda environment with software for running the pipeline

This installs [snakemake](https://snakemake.github.io/) and it's dependencies

```bash
mamba env create -f ./envs/pipeline_run_env.yaml
conda activate pipeline_run_env
```

## 10. Run the pipeline

First carry out a dry run

```bash
bash dryrun.sh
```

If there are no issues, start a full run

```bash
bash run.sh
```

## 11. Evaluate the pipeline run

Generate an interactive html report

```bash
bash report.sh
```

## 12. Commit and push to your forked version of the github repo

To maintain reproducibility, commit and push:

- All documentation
- All configuration files
- All run scripts
- The final report

## 13. Repeat step 12 each time you re-run the analysis with different parameters

## 14. Raise issues, create feature requests or create a pull request with the [upstream repo](https://github.com/ESR-NZ/human_genomics_pipeline) to merge any useful changes to the pipeline (optional)

See [the README](https://github.com/ESR-NZ/human_genomics_pipeline/blob/dev/README.md#contribute-back) for info on how to contribute back to the pipeline!
