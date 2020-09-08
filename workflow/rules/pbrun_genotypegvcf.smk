rule pbrun_genotypegvcf:
    input:
        gvcf = "../results/called/{family}_raw_snps_indels_tmp_combined.g.vcf",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output:
        protected("../results/called/{family}_raw_snps_indels.g.vcf")
    resources:
        gpu = 1
    log:
        "logs/pbrun_genotypegvcf/{family}.log"
    benchmark:
        "benchmarks/pbrun_genotypegvcf/{family}.tsv"
    message:
        "Performing joint genotyping on one or more samples pre-called with HaplotypeCaller for {input.gvcf}"
    shell:
        "pbrun genotypegvcf --in-gvcf {input.gvcf} --ref {input.refgenome} --out-vcf {output} &> {log}"
