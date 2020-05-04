rule multiqc_pre_trim:
    input:
        expand(["qc/fastqc/{sample}_R1_fastqc.zip", "qc/fastqc/{sample}_R2_fastqc.zip"], sample = SAMPLES)
    output:
        report("qc/multiqc/pre_trim_multiqc_report.html", caption = "../report/qualitychecks.rst", category = "Quality checks")
    conda:
        "../envs/multiqc.yaml"
    message:
        "Searching for analysis logs to compile a HTML report"
    shell:
        "multiqc {input} -o qc/multiqc/ -i pre_trim"