rule gatk_ApplyBQSR:
    input:
        bam = "../results/mapped/{sample}_sorted_mkdups.bam",
        recal_table = "../results/mapped/{sample}_recalibration_report.grp",
        ref = expand("{ref}", ref = config['REFGENOME']),
        dict = expand("{dict}", dict = config['DICT'])
    output:
        bam = protected("../results/mapped/{sample}_recalibrated.bam")
    params:
        java_opts = "-Xmx30g"
    log:
        "logs/gatk_ApplyBQSR/{sample}.log"
    benchmark:
        "benchmarks/gatk_ApplyBQSR/{sample}.tsv"
    conda:
        "../envs/gatk4.yaml"
    message:
        "Applying base quality score recalibration and producing a recalibrated BAM file for {input.bam}"
    wrapper:
        "0.64.0/bio/gatk/applybqsr"