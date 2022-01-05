# human_genomics_pipeline

A Snakemake workflow to process single samples (unrelated individuals) or cohorts (related individuals) of paired-end sequencing data (WGS or WES) using [bwa](http://bio-bwa.sourceforge.net/) and [GATK4](https://gatk.broadinstitute.org/hc/en-us). Quality control checks are also undertaken. The fastq files can be optionally trimmed with [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) and the pipeline can run on [NVIDIA GPU's](https://www.nvidia.com/en-gb/graphics-cards/) where [nvidia clara parabricks software is available](https://www.nvidia.com/en-us/docs/parabricks/quickstart-guide/software-overview/) for *significant* speedups in analysis times. This workflow is designed to follow the [GATK best practice workflow for germline short variant discovery (SNPs + Indels)](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-). This pipeline is designed to be followed by [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline) and the data ingested into [scout](https://github.com/Clinical-Genomics/scout) for clinical interpretation. However, this pipeline also stands on it's own, taking the data from fastq to vcf (raw sequencing data to called variants). This pipeline has been developed with human genetic data in mind, however we designed it to be species agnostic. Genetic data from other species can be analysed by setting a species-specific reference genome and variant databases in the configuration file (but not all situations have been tested).

- [human_genomics_pipeline](#human_genomics_pipeline)
  - [Workflow diagram - single samples](#workflow-diagram---single-samples)
  - [Workflow diagram - single samples - GPU accelerated](#workflow-diagram---single-samples---gpu-accelerated)
  - [Workflow diagram - cohort samples](#workflow-diagram---cohort-samples)
  - [Workflow diagram - cohort samples - GPU accelerated](#workflow-diagram---cohort-samples---gpu-accelerated)
  - [Prerequisites](#prerequisites)
  - [Test human_genomics_pipeline](#test-human_genomics_pipeline)
  - [Run human_genomics_pipeline](#run-human_genomics_pipeline)
  - [Contribute back!](#contribute-back)

## Workflow diagram - single samples

<img src="./images/rulegraph_single.png" class="center">

## Workflow diagram - single samples - GPU accelerated

<img src="./images/rulegraph_single_gpu.png" class="center">

## Workflow diagram - cohort samples

<img src="./images/rulegraph_cohort.png" class="center">

## Workflow diagram - cohort samples - GPU accelerated

<img src="./images/rulegraph_cohort_gpu.png" class="center">

## Prerequisites

- **Prerequisite hardware:** [NVIDIA GPUs](https://www.nvidia.com/en-gb/graphics-cards/) (for GPU accelerated runs)
- **Prerequisite software:** [NVIDIA CLARA PARABRICKS and dependencies](https://www.nvidia.com/en-us/docs/parabricks/local-installation/) (for GPU accelerated runs), [Git](https://git-scm.com/) (tested with version 2.7.4), [Mamba](https://github.com/TheSnakePit/mamba) (tested with version 0.4.4) with [Conda](https://docs.conda.io/projects/conda/en/latest/index.html) (tested with version 4.8.2), [gsutil](https://pypi.org/project/gsutil/) (tested with version 4.52), [gunzip](https://linux.die.net/man/1/gunzip) (tested with version 1.6)

## Test human_genomics_pipeline

The provided [test dataset](./test) can be used to test running this pipeline on a new machine, or test pipeline developments/releases. See instructions for using the test dataset on a [single machine like a laptop or single server/computer](./docs/running_on_a_single_machine.md#test-data) or a [high performance cluster](./docs/running_on_a_hpc.md#test-data).

## Run human_genomics_pipeline

See the docs for a walkthrough guide for running [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) on:

- [A single machine like a laptop or single server/computer](./docs/running_on_a_single_machine.md)
- [A high performance cluster](./docs/running_on_a_hpc.md)

## Contribute back!

- Raise issues in [the issues page](https://github.com/ESR-NZ/human_genomics_pipeline/issues)
- Create feature requests in [the issues page](https://github.com/ESR-NZ/human_genomics_pipeline/issues)
- Contribute your code! Create your own branch from the [development branch](https://github.com/ESR-NZ/human_genomics_pipeline/tree/dev) and create a pull request to the [development branch](https://github.com/ESR-NZ/human_genomics_pipeline/tree/dev) once the code is on point!

Contributions and feedback are always welcome! :blush:
