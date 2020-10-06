#!/bin/bash -x

snakemake \
--cores 32 \
--resources mem_mb=150000 \
--resources gpu=2 \
--use-conda \
--conda-frontend mamba \
--configfile ../config/config.yaml