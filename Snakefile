"""
Author: Miles Benton and Leah Kemp
Affiliation: ESR
Aim: A simple Snakemake workflow to process paired-end sequencing data (WGS) using bwa/GATK4. Designed to be used before vcf_annotation_pipeline.
Date created: 2019-08-21
Modified: 2020-04-29
Dry run: snakemake -n -j 24 --use-conda --configfile config.yaml
Full run: snakemake -j 24 --use-conda --configfile config.yaml
Report: snakemake --report report.html --configfile config.yaml --report-stylesheet stylesheet.css
Rule diagram: snakemake --rulegraph --configfile config.yaml | dot -Tpng > rulegraph.png
Workflow diagram (specific experiment): snakemake --dag --configfile config.yaml | dot -Tpng > dag.png
"""

##### Set up #####

# Define single or cohort samples from fastq dir using wildcards
if config['DATA'] == "Single":
    SAMPLES, = glob_wildcards("../fastq/{sample}_R1.fastq.gz")
elif config['DATA'] == "Cohort":
    # Fix cohort wildcard
    SAMPLES, = glob_wildcards("../fastq/{cohort}_R1.fastq.gz")
else: print("ERROR: Please check your configuration settings")

##### Target rules #####

rule all:
    input:
        expand("qc/multiqc/pre_trim_multiqc_report.html"),
        expand("qc/multiqc/post_trim_multiqc_report.html"),
        expand("vcf/{sample}_raw_snps_indels_AS_g.vcf", sample = SAMPLES)

#### Set up report #####

report: "report/workflow.rst"

##### load rules #####

localrules: multiqc_pre_trim, multiqc_post_trim, sambamba_index, sambamba_index_rgadd

include: "rules/fastqc.smk"
include: "rules/multiqc_pre_trim.smk"
include: "rules/trim_galore_pe.smk"
include: "rules/multiqc_post_trim.smk"
include: "rules/bwa_map.smk"
include: "rules/sambamba_sort.smk"
include: "rules/sambamba_mkdups.smk"
include: "rules/sambamba_index.smk"
include: "rules/gatk_add_replace_read_groups.smk"
include: "rules/sambamba_index_rgadd.smk"
include: "rules/gatk_base_recalibrator.smk"
include: "rules/gatk_apply_bqsr.smk"

if config['DATA'] == "Single":
    include: "rules/gatk_haplotype_caller_single.smk"

if config['DATA'] == "Cohort":
    include: "rules/gatk_haplotype_caller_gvcf.smk"
    include: "rules/gatk_combine_gvcf.smk"
    include: "rules/gatk_genotype_gvcf.smk"