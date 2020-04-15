rule multiqc_post_trim:
    input:
        zips = expand(["trim_galore/{sample}_R1_val_1_fastqc.zip", "trim_galore/{sample}_R2_val_2_fastqc.zip"], sample = SAMPLES)
    output:
        report("qc/post_trim_multiqc/post_trim_multiqc_report.html", caption = "../report/qualitychecks.rst", category = "Quality checks", subcategory = "Post-trimming")
    conda:
        "../envs/multiqc.yaml"
    message:
        "Searching for analysis logs to compile a HTML report"
    shell:
        "multiqc {input.zips} -o {output} -n post_trim_multiqc_report.html"