#!/bin/bash -x

snakemake \
--report ../results/report.html \
--configfile ../config/config.yaml \
--report-stylesheet ../config/ESR_stylesheet.css