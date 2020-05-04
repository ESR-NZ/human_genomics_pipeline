rule gatk4_recal:
    input:
        bams = "mapped/{sample}_sorted_mkdups_rgreplaced.bam",
        recal = "mapped/{sample}_recalibration_report.grp",
        genome = expand("{genome}", genome = config["FILEDIR"]["GENOME"])
    output:
        "mapped/{sample}_bwa_recal.bam"
    params:
        tdir = expand("{tdir}", tdir = config["TEMPDIR"])
    log:
        "logs/gatk_recal/{sample}.log"
    benchmark:
        "benchmarks/gatk_recal/{sample}.gatkrecal"
    conda:
        "../envs/gatk4.yaml"
    message:
	    "Applying base quality score recalibration and producing a recalibrated BAM file"
    shell:
        """
        gatk ApplyBQSR \
        -I {input.bams} \
        -bqsr {input.recal} \
        -R {input.genome} \
        -O {output} \
        --TMP_DIR {params.tdir}
        """