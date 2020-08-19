#!/bin/bash -x

snakemake \
-n \
-j 32 \
--use-conda \
--conda-frontend mamba \
--configfile ../config/config.yaml