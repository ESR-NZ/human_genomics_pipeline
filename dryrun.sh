#!/bin/bash -x

snakemake -n -j 24 --use-conda --configfile config.yaml