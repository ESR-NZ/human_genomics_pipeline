rule multiqc_pre_trim:
    input:
        dir = "qc/fastqc/"
    output:
        dir = "qc/pre_trim_multiqc/",
        rep = report("qc/multiqc_pre_trim/multiqc_report.html", caption = "../report/qualitychecks.rst", category = "Quality checks", subcategory = "Pre-trimming")
    conda:
        "../envs/multiqc.yaml"
    message:
        "Searching for analysis logs to compile a HTML report"
    shell:
        "multiqc {input.dir} -o {output.dir}"