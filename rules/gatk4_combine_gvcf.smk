rule gatk4_CombineGVCFs:
    input:
        # Fix wildcard naming of three input vcfs
        vcf1 = "vcf/{sample}_haplotype.vcf",
        vcf2 = "vcf/{sample}_haplotype.vcf",
        vcf3 = "vcf/{sample}_haplotype.vcf",
        genome = expand("{genome}", genome = config['FILEDIR']["GENOME"]),
        dbsnp = expand("{dbsnp}", dbsnp = config['FILEDIR']["dbSNP"])
    output:
        temp("vcf/{sample}_haplotype_gvcf_combined.vcf")
    params:
        "-G StandardAnnotation -G AS_StandardAnnotation"
    log:
        "logs/gatk_combinegvcf/{sample}.log"
    benchmark:
        "benchmarks/gatk_combinegvcf/{sample}.gatkcombinegvcf"
    conda:
        "../envs/gatk4.yaml"
    message:
        "Merging of GVCFs"
    shell: 
        """
        gatk --java-options "-Xmx64g -Xms64g" CombineGVCFs \
        -R {input.genome} --variant {input.vcf1} --variant {input.vcf2} --variant {input.vcf3} -O {output} {params}
        """