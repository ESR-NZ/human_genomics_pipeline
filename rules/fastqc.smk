rule fastqc:
    input:
        R1 = "../fastq/{sample}_R1.fastq.gz", sample = SAMPLES,
        R2 = "../fastq/{sample}_R2.fastq.gz", sample = SAMPLES
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
        "../envs/fastqc.yaml"
    shell:
        "fastqc {params} {input.R1} {input.R2} --outdir qc/fastqc 2> {log}"