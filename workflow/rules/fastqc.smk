rule fastqc:
    input:
        ["../../fastq/{sample}_1.fastq.gz", "../../fastq/{sample}_2.fastq.gz"]
    output:
        html = ["../results/qc/fastqc/{sample}_1_fastqc.html", "../results/qc/fastqc/{sample}_2_fastqc.html"],
        zip = ["../results/qc/fastqc/{sample}_1_fastqc.zip", "../results/qc/fastqc/{sample}_2_fastqc.zip"]
    log:
        "logs/fastqc/{sample}.log"
    benchmark:
        "benchmarks/fastqc/{sample}.tsv"
    singularity:
        "docker://biocontainers/fastqc:v0.11.9_cv8"
    threads: config['THREADS']
    resources:
        partition = config['PARTITION']['CPU']
    message:
        "Undertaking quality control checks on raw sequence data for {input}"
    shell:
        'fastqc {input} '
        '-o ../results/qc/fastqc/ '
        '-t {threads} '
        '&> {log}'