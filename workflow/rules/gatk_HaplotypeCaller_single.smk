rule gatk_HaplotypeCaller_single:
    input:
        bams = "../results/mapped/{sample}_recalibrated.bam",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME']),
        dbsnp = expand("{dbsnp}", dbsnp = config['dbSNP'])
    output:
        protected("../results/called/{sample}_raw_snps_indels.vcf")
    params:
        tdir = config['TEMPDIR'],
        padding = get_wes_padding_command,
        intervals = get_wes_intervals_command
    log:
        "logs/gatk_HaplotypeCaller_single/{sample}.log"
    benchmark:
        "benchmarks/gatk_HaplotypeCaller_single/{sample}.tsv"
    singularity:
        "docker://broadinstitute/gatk:4.2.6.1"
    threads: 1
    resources:
        mem_mb = config['MAXMEMORY'],
        partition = config['PARTITION']['CPU']
    message:
        "Calling germline SNPs and indels via local re-assembly of haplotypes for {input.bams}"
    shell:
        'gatk HaplotypeCaller '
        '--java-options "-Xmx{resources.mem_mb}m" '
        '-I {input.bams} '
        '-R {input.refgenome} '
        '-D {input.dbsnp} '
        '-O {output} '
        '--tmp-dir {params.tdir} {params.padding} {params.intervals} '
        '&> {log}'