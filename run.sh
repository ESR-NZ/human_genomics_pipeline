#!/bin/bash -x

snakemake \
-j 32 \
--use-conda \
--configfile config.yaml