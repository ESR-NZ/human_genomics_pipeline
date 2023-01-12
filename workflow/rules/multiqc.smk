rule multiqc:
    input:
        expand(["../results/qc/fastqc/{sample}_1_fastqc.zip", "../results/qc/fastqc/{sample}_2_fastqc.zip"], sample = SAMPLES)
    output:
        report("../results/qc/multiqc_report.html", caption = "../report/quality_checks.rst", category = "Quality checks")
    log:
        "logs/multiqc/multiqc.log"
    benchmark:
        "benchmarks/multiqc/multiqc.tsv"
    singularity:
        "docker://ewels/multiqc:v1.12"
    threads: 1
    resources:
        partition = config['PARTITION']['CPU']
    message:
        "Compiling a HTML report for quality control checks on raw sequence data"
    shell:
        'multiqc '
        '{input} '
        '-o ../results/qc/ '
        '&> {log}'