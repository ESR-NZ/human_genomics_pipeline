rule gatk_GenotypeGVCFs:
    input:
        gvcf = "../results/called/{family}_raw_snps_indels_tmp_combined.g.vcf",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output:
        protected("../results/called/{family}_raw_snps_indels.vcf")
    params:
        tdir = config['TEMPDIR'],
        padding = get_wes_padding_command,
        intervals = get_wes_intervals_command,
        other = "-G StandardAnnotation -G AS_StandardAnnotation"
    log:
        "logs/gatk_GenotypeGVCFs/{family}.log"
    benchmark:
        "benchmarks/gatk_GenotypeGVCFs/{family}.tsv"
    singularity:
        "docker://broadinstitute/gatk:4.2.6.1"
    threads: 1
    resources:
        mem_mb = config['MAXMEMORY'],
        partition = config['PARTITION']['CPU']
    message:
        "Performing joint genotyping on one or more samples pre-called with HaplotypeCaller for {input.gvcf}"
    shell:
        'gatk GenotypeGVCFs '
        '--java-options "-Xmx{resources.mem_mb}m" '
        '-R {input.refgenome} '
        '-V {input.gvcf} '
        '-O {output} {params.padding} {params.intervals} {params.other} '
        '&> {log}'