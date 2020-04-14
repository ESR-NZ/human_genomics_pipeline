"""
Author: Miles Benton
Affiliation: ESR
Aim: A simple Snakemake workflow to process paired-end sequencing data (WGS) using bwa/GATK4. Designed to be used before vcf_annotation_pipeline.
Date created: 2019-08-21
Modified: 2020-03-24
Run: snakemake -n -r -j 24 -p --use-conda
Rule diagram: snakemake --rulegraph | dot -Tpng > rulegraph.png
Workflow diagram (specific experiment): snakemake --dag | dot -Tpng > dag.png
"""

##### Set up #####

# define samples from fastq dir using wildcards
SAMPLES, = glob_wildcards("../fastq/{sample}_R1.fastq.gz")

##### Target rules #####

rule all:
    input:
        expand("qc/fastqc/{sample}_R1_fastqc.html", sample = SAMPLES),
        expand("vcf/{sample}.raw.snps.indels.AS.g.vcf", sample = SAMPLES)

##### load rules #####

include: "rules/fastqc.smk"
include: "rules/multiqc_pre_trim.smk"
include: "rules/trim_galore_pe.smk"
include: "rules/bwa_map.smk"
include: "rules/sambamba_sort.smk"
include: "rules/sambamba_mkdups.smk"
include: "rules/sambamba_index.smk"
include: "rules/gatk4_readgroup_add.smk"
include: "rules/sambamba_index_rgadd.smk"
include: "rules/gatk4_recal_report.smk"
include: "rules/gatk4_recal.smk"
include: "rules/gatk4_haplotypecaller.smk"