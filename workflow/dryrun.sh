#!/bin/bash -x

snakemake \
--dryrun \
--use-conda \
--conda-frontend mamba \
--use-singularity \
--singularity-args '-B /bind/location/' \
--configfile ../config/config.yaml