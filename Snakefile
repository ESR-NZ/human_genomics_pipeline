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
# adapt paths as appropriate
# GRCh37
GENOME = "/store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta"
# GRCh38
# GENOME = "/data/publicData/genomes/human/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
# temp dir
TEMPDIR = "/store/lkemp/tmp/"
# dbSNP
dbSNP = "/store/lkemp/publicData/dbSNP/ncbi/GRCh37/build151/GATK/All_20180423.vcf.gz"

# define samples from fastq dir using wildcards
SAMPLES, = glob_wildcards("../data/exomes/fastq/{sample}_R1.fastq.gz")

rule all:
    input:
        expand("qc/fastqc/{sample}_R1_fastqc.html", sample=SAMPLES),
        "qc/pre_trim_multiqc/multiqc_report.html",
        expand("vcf/{sample}.raw.snps.indels.AS.g.vcf", sample = SAMPLES)

rule fastqc:
    input:
        R1 = "../data/exomes/fastq/{sample}_R1.fastq.gz",
        R2 = "../data/exomes/fastq/{sample}_R2.fastq.gz"
    output:
        html = ["qc/fastqc/{sample}_R1_fastqc.html", "qc/fastqc/{sample}_R2_fastqc.html"],
        zip = ["qc/fastqc/{sample}_R1_fastqc.zip", "qc/fastqc/{sample}_R2_fastqc.zip"]
    params: 
        "--threads 4"
    log:
        "logs/fastqc/{sample}.log"
    benchmark:
        "benchmarks/fastqc/{sample}.fastqc"
    conda:
        "envs/fastqc.yaml"
    shell:
        "fastqc {params} {input.R1} {input.R2} --outdir qc/fastqc 2> {log}"

rule multiqc_pre_trim:
    input:
        zips = expand(["qc/fastqc/{sample}_R1_fastqc.zip", "qc/fastqc/{sample}_R2_fastqc.zip"], sample = SAMPLES)
    output:
        "qc/pre_trim_multiqc/"
    conda:
        "envs/multiqc.yaml"
    shell:
        "multiqc {input.zips} --outdir {output}"

rule trim_galore_pe:
    input:
        R1 = "../data/exomes/fastq/{sample}_R1.fastq.gz",
        R2 = "../data/exomes/fastq/{sample}_R2.fastq.gz"
    output:
        "trim_galore/{sample}_R1_val_1.fq.gz",
        "trim_galore/{sample}_R2_val_2.fq.gz"
    params:
        extra="--illumina -q 20"
    log:
        "logs/trim_galore/{sample}.log"
    benchmark:
        "benchmarks/trim_galore/{sample}.trim"
    conda:
        "envs/trim_galore.yaml"
    shell:
        "trim_galore --illumina --fastqc --paired {input.R1} {input.R2} --output_dir trim_galore/"

rule bwa_map:
    input:
        R1 = "trim_galore/{sample}_R1_val_1.fq.gz",
        R2 = "trim_galore/{sample}_R2_val_2.fq.gz"
    output: 
        temp("mapped/{sample}_bwamem.bam")
    log:
        "logs/bwamem/{sample}.log"
    benchmark:
        "benchmarks/bwamem/{sample}.bwamem"
    conda:
        "envs/bwa.yaml"
    # params:
    #     ""
    threads: 12
    shell: "bwa mem -M -t {threads} {GENOME} {input.R1} {input.R2} | samtools view -@ {threads} -Sbh - > {output}"

rule sambamba_sort:
    input:
        bams = "mapped/{sample}_bwamem.bam"
    output:
        temp("mapped/{sample}_bwamem_sorted.bam")
    log:
        "logs/sambamba_sort/{sample}.log"
    benchmark:
        "benchmarks/sambamba_sort/{sample}.sambamba"
    conda:
        "envs/sambamba.yaml"
    threads: 4
    shell:
        "sambamba sort -t {threads} -m 6G --tmpdir={TEMPDIR} -p -o {output} {input.bams}"

rule sambamba_mkdups:
    input:
        bams = "mapped/{sample}_bwamem_sorted.bam"
    output:
        temp("mapped/{sample}_bwamem_sorted_mkdups.bam")
    log:
        "logs/sambamba_mkdups/{sample}.log"
    benchmark:
        "benchmarks/sambamba_mkdups/{sample}.sambamba"
    conda:
        "envs/sambamba.yaml"
    threads: 4
    params: "--sort-buffer-size=6144 --overflow-list-size=600000 --hash-table-size=600000"
    shell:
        "sambamba markdup -t {threads} {params} --tmpdir={TEMPDIR} -p {input.bams} {output}"

rule sambamba_index:
    input:
        bams = "mapped/{sample}_bwamem_sorted_mkdups.bam"
    output:
        temp("mapped/{sample}_bwamem_sorted_mkdups.bam.bai")
    log:
        "logs/sambamba_index/{sample}.log"
    benchmark:
        "benchmarks/sambamba_index/{sample}.sambamba"
    conda:
        "envs/sambamba.yaml"
    threads: 4
    shell:
        "sambamba index -t {threads} -p {input.bams}"

rule gatk4_readgroup_add:
    input:
        bams = "mapped/{sample}_bwamem_sorted_mkdups.bam"
    output:
        temp("mapped/{sample}_sorted_mkdups_rgreplaced.bam")
    log:
        "logs/gatk_readgroup/{sample}.log"
    benchmark:
        "benchmarks/gatk_readgroup/{sample}.readgroup"
    conda:
        "envs/gatk4.yaml"
    threads: 4
    shell:
        "gatk AddOrReplaceReadGroups --INPUT {input.bams} --OUTPUT {output} --RGID 4 --RGLB lib1 --RGPL illumina --RGPU unit1 --RGSM 20"

rule sambamba_index_rgadd:
    input:
        bams = "mapped/{sample}_sorted_mkdups_rgreplaced.bam"
    output:
        temp("mapped/{sample}_sorted_mkdups_rgreplaced.bam.bai")
    log:
        "logs/sambamba_index/{sample}_rg.log"
    benchmark:
        "benchmarks/sambamba_index/{sample}_rg.sambamba"
    conda:
        "envs/sambamba.yaml"
    threads: 4
    shell:
        "sambamba index -t {threads} -p {input.bams}"

rule gatk4_recal_report:
    input:
        bams = "mapped/{sample}_sorted_mkdups_rgreplaced.bam"
    output:
        "mapped/{sample}_recalibration_report.grp"
    log:
        "logs/gatk_recalrep/{sample}.log"
    benchmark:
        "benchmarks/gatk_recalrep/{sample}.gatkrecalrep"
    conda:
        "envs/gatk4.yaml"
    threads: 4
    shell:
        "gatk BaseRecalibrator --reference {GENOME} --input {input.bams} --known-sites {dbSNP} --output {output}"

rule gatk4_recal:
    input:
        bams = "mapped/{sample}_sorted_mkdups_rgreplaced.bam",
        recal = "mapped/{sample}_recalibration_report.grp"
    output:
        "mapped/{sample}_bwa_recal.bam"
    log:
        "logs/gatk_recal/{sample}.log"
    benchmark:
        "benchmarks/gatk_recal/{sample}.gatkrecal"
    conda:
        "envs/gatk4.yaml"
    threads: 4
    shell:
        "gatk ApplyBQSR --reference {GENOME} --bqsr-recal-file {input.recal} --input {input.bams} --output {output}"

rule gatk4_HaplotypeCaller:
    input:
        bams = "mapped/{sample}_bwa_recal.bam"
    output:
        "vcf/{sample}.raw.snps.indels.AS.g.vcf"
    log:
        "logs/gatk_haplocall/{sample}.log"
    benchmark:
        "benchmarks/gatk_haplocall/{sample}.gatkhaplocall"
    conda:
        "envs/gatk4.yaml"
    threads: 4
    shell:
        "gatk HaplotypeCaller --reference {GENOME} --emit-ref-confidence GVCF --dbsnp {dbSNP} --input {input.bams} --output {output} --tmp-dir {TEMPDIR}"