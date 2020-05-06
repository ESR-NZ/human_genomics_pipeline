rule gatk4_CombineGVCFs:
    input:
        # Fix wildcard naming of three input vcfs
        vcf1 = "vcf/{cohort}_{sample}_haplotype.vcf",
        vcf2 = "vcf/{cohort}_{sample}_haplotype.vcf",
        vcf3 = "vcf/{cohort}_{sample}_haplotype.vcf",
        refgenome = expand("{refgenome}", refgenome = config['GENOME']),
        dbsnp = expand("{dbsnp}", dbsnp = config['dbSNP'])
    output:
        temp("vcf/{sample}_haplotype_gvcf_combined.vcf")
    params:
        padding = expand("{padding}", padding = config['WES']['PADDING']),
        intervals = expand("{intervals}", intervals = config['WES']['INTERVALS']),
        other = "-G StandardAnnotation -G AS_StandardAnnotation"
    log:
        "logs/gatk_combinegvcf/{sample}.log"
    benchmark:
        report("benchmarks/gatk_combinegvcf/{sample}.gatkcombinegvcf", caption = "../report/benchmarking.rst", category = "Benchmarking")
    conda:
        "../envs/gatk4.yaml"
    message:
        "Merging of GVCFs"
    shell: 
        """
        gatk --java-options "-Xmx64g -Xms64g" CombineGVCFs \
        -R {input.refgenome} --variant {input.vcf1} --variant {input.vcf2} --variant {input.vcf3} -O {output} {params.padding} {params.intervals} {params.other}
        """