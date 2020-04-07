rule multiqc_pre_trim:
    input:
        zips = expand(["qc/fastqc/{sample}_R1_fastqc.zip", "qc/fastqc/{sample}_R2_fastqc.zip"], sample = SAMPLES)
    output:
        "qc/pre_trim_multiqc/"
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc {input.zips} --outdir {output}"