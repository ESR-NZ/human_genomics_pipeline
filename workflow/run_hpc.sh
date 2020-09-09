#!/bin/bash -x

snakemake \
-j 32 \
--resources gpu=2 \
--use-conda \
--conda-frontend mamba \
--configfile ../config/config.yaml \
--cluster-config ../config/cluster.json \
--cluster "sbatch -A {cluster.account} \
-p {cluster.partition} \
--nodes {cluster.nodes} \
--ntasks {cluster.ntasks}"