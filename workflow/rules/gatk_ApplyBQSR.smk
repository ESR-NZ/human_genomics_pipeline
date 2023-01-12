rule gatk_ApplyBQSR:
    input:
        bam = "../results/mapped/{sample}_sorted_mkdups.bam",
        recal = "../results/mapped/{sample}_recalibration_report.grp",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output:
        bam = protected("../results/mapped/{sample}_recalibrated.bam")
    params:
        padding = get_wes_padding_command,
        intervals = get_wes_intervals_command
    log:
        "logs/gatk_ApplyBQSR/{sample}.log"
    benchmark:
        "benchmarks/gatk_ApplyBQSR/{sample}.tsv"
    singularity:
        "docker://broadinstitute/gatk:4.2.6.1"
    threads: 1
    resources:
        mem_mb = config['MAXMEMORY'],
        partition = config['PARTITION']['CPU']
    message:
        "Applying base quality score recalibration and producing a recalibrated BAM file for {input.bam}"
    shell:
        'gatk ApplyBQSR '
        '--java-options "-Xmx{resources.mem_mb}m" '
        '-I {input.bam} '
        '-bqsr {input.recal} '
        '-R {input.refgenome} '
        '-O {output} {params.padding} {params.intervals}'