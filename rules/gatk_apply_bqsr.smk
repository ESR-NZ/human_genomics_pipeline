rule gatk4_ApplyBQSR:
    input:
        bams = "mapped/{sample}_sorted_mkdups_rgreplaced.bam",
        index = "mapped/{sample}_sorted_mkdups_rgreplaced.bam.bai",
        recal = "mapped/{sample}_recalibration_report.grp",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output:
        "mapped/{sample}_bwa_recal.bam"
    params:
        padding = expand("{padding}", padding = config['WES']['PADDING']),
        intervals = expand("{intervals}", intervals = config['WES']['INTERVALS'])
    log:
        "logs/gatk_recal/{sample}.log"
    benchmark:
        "benchmarks/gatk_recal/{sample}.gatkrecal"
    conda:
        "../envs/gatk4.yaml"
    message:
        "Applying base quality score recalibration and producing a recalibrated BAM file"
    shell:
        "gatk ApplyBQSR -I {input.bams} -bqsr {input.recal} -R {input.refgenome} -O {output} {params.padding} {params.intervals}"