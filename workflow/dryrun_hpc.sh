#!/bin/bash -x

snakemake \
--dryrun \
--profile ../config/slurm \
# --singularity-args '-B /bind/location/' \
