rule pbrun_triocombinegvcf:
    input:
        vcf_dummy = expand("../results/called/{sample}_raw_snps_indels_tmp.g.vcf", sample = SAMPLES), # a dummy vcf to connect this rule to gatk_HaplotypeCaller
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output:
        temp("../results/called/{family}_raw_snps_indels_tmp_combined.g.vcf")
    params:
        command = get_parabricks_combinegvcf_command
    log:
        "logs/pbrun_triocombinegvcf/{family}.log"
    benchmark:
        "benchmarks/pbrun_triocombinegvcf/{family}.tsv"
    threads: config['THREADS']
    resources:
        gpu = config['GPU'],
        partition = config['PARTITION']['GPU'],
        job_name = "pbrun_triocombinegvcf_gpu"
    message:
        "Merging one or more HaplotypeCaller GVCF files into a single GVCF"
    shell:
        'pbrun triocombinegvcf '
        '--ref {input.refgenome} '
        '{params.command} '
        '--out-variants {output} '
        '&> {log}'