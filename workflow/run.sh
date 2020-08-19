#!/bin/bash -x

snakemake \
-j 32 \
--use-conda \
--conda-frontend mamba \
--configfile ../config/config.yaml