rule multiqc_pre_trim:
    input:
        dir = "qc/fastqc/"
    output:
        dir = "qc/pre_trim_multiqc/"
    conda:
        "../envs/multiqc.yaml"
    message:
        "Searching for analysis logs to compile a HTML report"
    shell:
        "multiqc {input.dir} -o {output.dir}"