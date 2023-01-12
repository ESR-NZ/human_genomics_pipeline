#!/bin/bash -x
snakemake \
--profile ../config/slurm/ \
# --singularity-args '-B /bind/location/' \ #TODO: Find out what this is
--configfile ../config/config.yaml