#!/bin/bash -x

snakemake \
-j 24 \
--use-conda \
--configfile config.yaml