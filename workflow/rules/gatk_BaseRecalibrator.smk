rule gatk_BaseRecalibrator:
    input:
        bams = "../results/mapped/{sample}_sorted_mkdups.bam",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output:
        report("../results/mapped/{sample}_recalibration_report.grp", caption = "../report/recalibration.rst", category = "Base recalibration")
    params:
        tdir = expand("{tdir}", tdir = config['TEMPDIR']),
        padding = expand("{padding}", padding = config['WES']['PADDING']),
        intervals = expand("{intervals}", intervals = config['WES']['INTERVALS']),
        recalibration_resources = expand("{recalibration_resources}", recalibration_resources = config['RECALIBRATION']['RESOURCES']),
        other = "--java-options -Xmx30g"
    log:
        "logs/gatk_BaseRecalibrator/{sample}.log"
    benchmark:
        "benchmarks/gatk_BaseRecalibrator/{sample}.tsv"
    conda:
        "../envs/gatk4.yaml"
    message:
        "Generating a recalibration table for {input.bams}"
    shell:
        "gatk BaseRecalibrator {params.other} -I {input.bams} -R {input.refgenome} -O {output} --tmp-dir {params.tdir} {params.padding} {params.intervals} {params.recalibration_resources} &> {log}"