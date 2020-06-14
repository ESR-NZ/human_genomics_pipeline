rule trim_galore_pe:
    input:
        R1 = "../fastq/{sample}_R1.fastq.gz",
        R2 = "../fastq/{sample}_R2.fastq.gz"
    output:
        "trim_galore/{sample}_R1_val_1.fq.gz",
        "trim_galore/{sample}_R2_val_2.fq.gz",
        "trim_galore/{sample}_R1_val_1_fastqc.zip",
        "trim_galore/{sample}_R2_val_2_fastqc.zip",
        report("trim_galore/{sample}_R1.fastq.gz_trimming_report.txt", caption = "../report/trimming_R1.rst", category = "Trimming"),
        report("trim_galore/{sample}_R2.fastq.gz_trimming_report.txt", caption = "../report/trimming_R2.rst", category = "Trimming")
    params:
        "--illumina --fastqc -q 20"
    log:
        "logs/trim_galore/{sample}.log"
    benchmark:
        "benchmarks/trim_galore/{sample}.trim"
    conda:
        "../envs/trim_galore.yaml"
    threads: 16
    message:
        "Applying quality and adapter trimming of input fastq files: {input.R1} and {input.R2}"
    shell:
        "trim_galore --paired {input.R1} {input.R2} -o trim_galore/ {params} --cores {threads}"