rule multiqc:
    input:
        fastqc = expand("../results/qc/fastqc/{sample}_fastqc.zip", sample = SAMPLES),
        trimming1 = expand("../results/trimmed/{sample}_1.fastq.gz_trimming_report.txt", sample = SAMPLES),
        trimming2 = expand("../results/trimmed/{sample}_2.fastq.gz_trimming_report.txt", sample = SAMPLES)
    output:
        report("../results/qc/multiqc_report.html", caption = "../report/quality_checks.rst", category = "Quality checks")
    conda:
        "../envs/multiqc.yaml"
    log:
        "logs/multiqc/multiqc.log"
    message:
        "Compiling a HTML report for quality control checks on raw sequence data"
    shell:
        "multiqc {input.fastqc} {input.trimming1} {input.trimming2} -o ../results/qc/ &> {log}"
