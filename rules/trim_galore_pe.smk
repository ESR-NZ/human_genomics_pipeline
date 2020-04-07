rule trim_galore_pe:
    input:
        R1 = "../fastq/{sample}_R1.fastq.gz", sample = SAMPLES,
        R2 = "../fastq/{sample}_R2.fastq.gz", sample = SAMPLES
    output:
        "trim_galore/{sample}_R1_val_1.fq.gz",
        "trim_galore/{sample}_R2_val_2.fq.gz"
    params:
        extra = "--illumina -q 20"
    log:
        "logs/trim_galore/{sample}.log"
    benchmark:
        "benchmarks/trim_galore/{sample}.trim"
    conda:
        "../envs/trim_galore.yaml"
    shell:
        "trim_galore --illumina --fastqc --paired {input.R1} {input.R2} --output_dir trim_galore/"