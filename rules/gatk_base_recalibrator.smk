rule gatk4_BaseRecalibrator:
    input:
        bams = "mapped/{sample}_sorted_mkdups_rgreplaced.bam",
        index = "mapped/{sample}_sorted_mkdups_rgreplaced.bam.bai",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output:
        report("mapped/{sample}_recalibration_report.grp", caption = "../report/recalibration.rst", category = "Base recalibration")
    params:
        padding = expand("{padding}", padding = config['WES']['PADDING']),
        intervals = expand("{intervals}", intervals = config['WES']['INTERVALS']),
        recalibration_resources = expand("{recalibration_resources}", recalibration_resources = config['RECALIBRATION']['RESOURCES'])
    log:
        "logs/gatk_recalrep/{sample}.log"
    benchmark:
        "benchmarks/gatk_recalrep/{sample}.gatkrecalrep"
    conda:
        "../envs/gatk4.yaml"
    message:
        "Generating a recalibration table for the following rule (Base Quality Score Recalibration)"
    shell:
        "gatk BaseRecalibrator -I {input.bams} -R {input.refgenome} -O {output} {params.padding} {params.intervals} {params.recalibration_resources}"