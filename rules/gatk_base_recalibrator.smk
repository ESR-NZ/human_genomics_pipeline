rule gatk4_BaseRecalibrator:
    input:
        bams = "mapped/{sample}_sorted_mkdups_rgreplaced.bam",
        index = "mapped/{sample}_sorted_mkdups_rgreplaced.bam.bai",
        genome = expand("{genome}", genome = config['FILEDIR']['GENOME']),
        dbsnp = expand("{dbsnp}", dbsnp = config['FILEDIR']['dbSNP'])
    output:
        report("mapped/{sample}_recalibration_report.grp", caption = "../report/recalibration.rst", category = "Base recalibration")
    params:
        padding = expand("{padding}", padding = config['WES']['PADDING']),
        intervals = expand("{intervals}", intervals = config['WES']['INTERVALS'])
    log:
        "logs/gatk_recalrep/{sample}.log"
    benchmark:
        report("benchmarks/gatk_recalrep/{sample}.gatkrecalrep", caption = "../report/benchmarking.rst", category = "Benchmarking")
    conda:
        "../envs/gatk4.yaml"
    message:
        "Generating a recalibration table for the following rule (Base Quality Score Recalibration)"
    shell:
        "gatk BaseRecalibrator -I {input.bams} -R {input.genome} --known-sites {input.dbsnp} -O {output} {params.padding} {params.intervals}"