#!/bin/bash -x

snakemake \
--use-conda \
--conda-frontend mamba \
--latency-wait 120 \
--use-singularity \
--profile ../config/slurm/ \
--singularity-args '-B /bind/location/' \
--configfile ../config/config.yaml