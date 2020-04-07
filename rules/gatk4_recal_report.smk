rule gatk4_recal_report:
    input:
        bams = "mapped/{sample}_sorted_mkdups_rgreplaced.bam"
    output:
        "mapped/{sample}_recalibration_report.grp"
    params:
        genome = expand("{genome}", genome = config["GENOME"]),
        dbsnp = expand("{dbsnp}", dbsnp = config["dbSNP"])
    log:
        "logs/gatk_recalrep/{sample}.log"
    benchmark:
        "benchmarks/gatk_recalrep/{sample}.gatkrecalrep"
    conda:
        "../envs/gatk4.yaml"
    threads: 4
    message:
	 "Generating a recalibration table for the following rule (Base Quality Score Recalibration)"
    shell:
        "gatk BaseRecalibrator --reference {params.genome} --input {input.bams} --known-sites {params.dbsnp} --output {output}"