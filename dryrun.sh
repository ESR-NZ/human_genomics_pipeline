#!/bin/bash -x

snakemake \
-n -j 32 \
--use-conda \
--configfile config.yaml