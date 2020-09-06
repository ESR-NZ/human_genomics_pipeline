rule multiqc:
    input:
        expand("../results/qc/fastqc/{sample}_fastqc.zip", sample = SAMPLES)
    output:
        report("../results/qc/multiqc_report.html", caption = "../report/quality_checks.rst", category = "Quality checks")
    conda:
        "../envs/multiqc.yaml"
    log:
        "logs/multiqc/multiqc.log"
    message:
        "Compiling a HTML report for quality control checks on raw sequence data"
    shell:
        "multiqc {input} -o ../results/qc/ &> {log}"