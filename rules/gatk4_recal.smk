rule gatk4_recal:
    input:
        bams = "mapped/{sample}_sorted_mkdups_rgreplaced.bam",
        recal = "mapped/{sample}_recalibration_report.grp"
    output:
        "mapped/{sample}_bwa_recal.bam"
    params:
        genome = expand("{genome}", genome = config["GENOME"])
    log:
        "logs/gatk_recal/{sample}.log"
    benchmark:
        "benchmarks/gatk_recal/{sample}.gatkrecal"
    conda:
        "../envs/gatk4.yaml"
    message:
	    "Applying base quality score recalibration and producing a recalibrated BAM file"
    shell:
        "gatk ApplyBQSR -I {input.bams} -bqsr {input.recal} -O {output} -R {params.genome}"