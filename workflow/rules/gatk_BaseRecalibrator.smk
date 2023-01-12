rule gatk_BaseRecalibrator:
    input:
        bams = "../results/mapped/{sample}_sorted_mkdups.bam",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output:
        report("../results/mapped/{sample}_recalibration_report.grp", caption = "../report/recalibration.rst", category = "Base recalibration")
    params:
        tdir = config['TEMPDIR'],
        padding = get_wes_padding_command,
        intervals = get_wes_intervals_command,
        recalibration_resources = get_recal_resources_command
    log:
        "logs/gatk_BaseRecalibrator/{sample}.log"
    benchmark:
        "benchmarks/gatk_BaseRecalibrator/{sample}.tsv"
    singularity:
        "docker://broadinstitute/gatk:4.2.6.1"
    threads: 1
    resources:
        mem_mb = config['MAXMEMORY'],
        partition = config['PARTITION']['CPU']
    message:
        "Generating a recalibration table for {input.bams}"
    shell:
        'gatk BaseRecalibrator '
        '--java-options "-Xmx{resources.mem_mb}m" '
        '-I {input.bams} '
        '-R {input.refgenome} '
        '-O {output} '
        '--tmp-dir {params.tdir} {params.padding} {params.intervals} {params.recalibration_resources} '
        '&> {log}'