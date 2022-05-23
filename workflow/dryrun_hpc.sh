#!/bin/bash -x

snakemake \
--dryrun \
--use-conda \
--conda-frontend mamba \
--latency-wait 120 \
--use-singularity \
--profile ../config/slurm/ \
--singularity-args '-B /bind/location/' \
--configfile ../config/config.yaml