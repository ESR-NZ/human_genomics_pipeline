rule gatk4_BaseRecalibrator:
    input:
        bams = "mapped/{sample}_sorted_mkdups_rgreplaced.bam",
        genome = expand("{genome}", genome = config['FILEDIR']['GENOME']),
        dbsnp = expand("{dbsnp}", dbsnp = config['FILEDIR']['dbSNP'])
    output:
        report("mapped/{sample}_recalibration_report.grp", caption = "../report/recalibration.rst", category = "Base recalibration")
    log:
        "logs/gatk_recalrep/{sample}.log"
    benchmark:
        "benchmarks/gatk_recalrep/{sample}.gatkrecalrep"
    conda:
        "../envs/gatk4.yaml"
    message:
        "Generating a recalibration table for the following rule (Base Quality Score Recalibration)"
    shell:
        "gatk BaseRecalibrator -I {input.bams} -R {input.genome} --known-sites {input.dbsnp} -O {output}"