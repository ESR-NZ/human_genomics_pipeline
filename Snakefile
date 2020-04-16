"""
Author: Miles Benton
Affiliation: ESR
Aim: A simple Snakemake workflow to process paired-end sequencing data (WGS) using bwa/GATK4. Designed to be used before vcf_annotation_pipeline.
Date created: 2019-08-21
Modified: 2020-04-14
Run: snakemake -n -r -j 24 -p --use-conda --configfile config_GRCh38.yaml
Rule diagram: snakemake --rulegraph --configfile config_GRCh38.yaml | dot -Tpng > rulegraph.png
Workflow diagram (specific experiment): snakemake --dag --configfile config_GRCh38.yaml | dot -Tpng > dag.png
"""

##### Set up #####

# define samples from fastq dir using wildcards
SAMPLES, = glob_wildcards("../fastq/{sample}_R1.fastq.gz")

##### Target rules #####

rule all:
    input:
        expand("qc/multiqc/pre_trim_multiqc_report.html"),
        expand("qc/multiqc/post_trim_multiqc_report.html"),
        expand("vcf/{sample}.raw.snps.indels.AS.g.vcf", sample = SAMPLES)

#### Set up report #####

report: config["REPORTWORKFLOW"]

##### load rules #####

include: "rules/fastqc.smk"
include: "rules/multiqc_pre_trim.smk"
include: "rules/trim_galore_pe.smk"
include: "rules/multiqc_post_trim.smk"
include: "rules/bwa_map.smk"
include: "rules/sambamba_sort.smk"
include: "rules/sambamba_mkdups.smk"
include: "rules/sambamba_index.smk"
include: "rules/gatk4_readgroup_add.smk"
include: "rules/sambamba_index_rgadd.smk"
include: "rules/gatk4_recal_report.smk"
include: "rules/gatk4_recal.smk"
include: "rules/gatk4_haplotypecaller.smk"