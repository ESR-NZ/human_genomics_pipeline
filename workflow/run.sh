#!/bin/bash -x

snakemake \
--latency-wait 600 \
--jobs 4 \
--resources mem_mb=8000 \ # Will conflict with the config.yaml
--threads 8 \ # Will conflict with the config.yaml
--use-conda \
--conda-frontend mamba \
--use-singularity \
# --singularity-args '-B /bind/location/' \
--configfile ../config/config.yaml
--keep-going \
--keep-incomplete \
--reason \
--printshellcmds \
--rerun-incomplete