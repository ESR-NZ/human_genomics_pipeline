rule fastqc:
    input:
        R1 = "../fastq/{sample}_R1.fastq.gz",
        R2 = "../fastq/{sample}_R2.fastq.gz"
    output:
        html = ["qc/fastqc/{sample}_R1_fastqc.html", "qc/fastqc/{sample}_R2_fastqc.html"],
        zip = ["qc/fastqc/{sample}_R1_fastqc.zip", "qc/fastqc/{sample}_R2_fastqc.zip"]
    log:
        "logs/fastqc/{sample}.log"
    benchmark:
        "benchmarks/fastqc/{sample}.fastqc"
    conda:
        "../envs/fastqc.yaml"
    threads: 2
    message:
        "Undertaking quality control checks on raw sequence data"
    shell:
        "fastqc {input.R1} {input.R2} --outdir qc/fastqc 2> {log} --threads {threads}"