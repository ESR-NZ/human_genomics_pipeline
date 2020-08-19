rule multiqc:
    input:
        expand("../results/qc/fastqc/{sample}_fastqc.zip", sample = SAMPLES)
    output:
        report("../results/qc/multiqc.html", caption = "../report/quality_checks.rst", category = "Quality checks")
    params:
        ""
    log:
        "logs/multiqc/multiqc.log"
    message:
        "Compiling a HTML report for quality control checks on raw sequence data"
    wrapper:
        "0.64.0/bio/multiqc"