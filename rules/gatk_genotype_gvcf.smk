rule gatk4_GenotypeGVCFs:
    input:
        vcf = "raw/{sample}_haplotype_combined.vcf",
        genome = expand("{genome}", genome = config['FILEDIR']['GENOME'])
    output:
        "vcf/{sample}_raw_snps_indels_AS_g.vcf"
    params:
        "-G StandardAnnotation -G AS_StandardAnnotation"
    log:
        "logs/gatk_genotype_gvcf/{sample}.log"
    benchmark:
        "benchmarks/gatk_genotype_gvcf/{sample}.gatkgenotypegvcf"
    conda:
        "../envs/gatk4.yaml"
    message:
        "Performing joint genotyping"
    shell:
        """
        gatk --java-options "-Xmx64g -Xms64g" GenotypeGVCFs \
            -R {input.genome} -V {input.vcf} -O {output}
        """