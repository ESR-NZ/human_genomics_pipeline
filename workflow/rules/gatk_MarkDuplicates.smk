rule gatk_MarkDuplicates:
    input:
        "../results/mapped/{sample}_sorted.bam"
    output:
        bam = temp("../results/mapped/{sample}_sorted_mkdups.bam"),
        metrics = "../results/mapped/{sample}_sorted_mkdups_metrics.txt"
    params:
        tdir = config['TEMPDIR']
    log:
        "logs/gatk_MarkDuplicates/{sample}.log"
    benchmark:
        "benchmarks/gatk_MarkDuplicates/{sample}.tsv"
    singularity:
        "docker://broadinstitute/gatk:4.2.6.1"
    threads: 1
    resources:
        mem_mb = config['MAXMEMORY'],
        partition = config['PARTITION']['CPU']
    message:
        "Locating and tagging duplicate reads in {input}"
    shell:
        'gatk MarkDuplicates '
        '--java-options "-Xmx{resources.mem_mb}m" '
        '-I {input} '
        '-O {output.bam} '
        '-M {output.metrics} '
        '--TMP_DIR {params.tdir} '
        '&> {log}'