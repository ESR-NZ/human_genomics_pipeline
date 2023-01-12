rule gatk_CombineGVCFs:
    input:
        vcf_dummy = expand("../results/called/{sample}_raw_snps_indels_tmp.g.vcf", sample = SAMPLES), # a dummy vcf to connect this rule to gatk_HaplotypeCaller
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output: 
        vcf = temp("../results/called/{family}_raw_snps_indels_tmp_combined.g.vcf"),
        index = temp("../results/called/{family}_raw_snps_indels_tmp_combined.g.vcf.idx")
    params:
        command = get_gatk_combinegvcf_command,
        tdir = config['TEMPDIR'],
        padding = get_wes_padding_command,
        intervals = get_wes_intervals_command,
        other = "-G StandardAnnotation -G AS_StandardAnnotation"
    log: 
        "logs/gatk_CombineGVCFs/{family}.log"
    benchmark:
        "benchmarks/gatk_CombineGVCFs/{family}.tsv"
    singularity:
        "docker://broadinstitute/gatk:4.2.6.1"
    threads: 1
    resources:
        partition = config['PARTITION']['CPU']
    message:
        "Merging one or more HaplotypeCaller GVCF files into a single GVCF"
    shell:
        'gatk CombineGVCFs '
        '-R {input.refgenome} {params.command} '
        '-O {output.vcf} 
        '--tmp-dir {params.tdir} {params.other} '
        '&> {log}'