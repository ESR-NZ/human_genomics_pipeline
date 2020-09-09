#!/bin/bash -x

snakemake \
-n \
-j 32 \
--resources gpu=2 \
--use-conda \
--conda-frontend mamba \
--configfile ../config/config.yaml