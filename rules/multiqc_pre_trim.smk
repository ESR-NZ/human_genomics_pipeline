rule multiqc_pre_trim:
    input:
        zips = expand(["qc/fastqc/{sample}_R1_fastqc.zip", "qc/fastqc/{sample}_R2_fastqc.zip"], sample = SAMPLES)
    output:
        dir = "qc/pre_trim_multiqc/",
        rep = report("qc/pre_trim_multiqc/multiqc_report.html", caption = "../report/qualitychecks.rst", category = "Quality checks", subcategory = "Pre-trimming")
    conda:
        "../envs/multiqc.yaml"
    message:
        "Searching for analysis logs to compile a HTML report"
    shell:
        "multiqc {input.zips} -o {output.dir}"