rule gatk4_GenotypeGVCFs:
    input:
        vcf = "vcf/{sample}_haplotype_gvcf_combined.vcf",
        genome = expand("{genome}", genome = config['FILEDIR']['GENOME'])
    output:
        "vcf/{sample}_raw_snps_indels_AS_g.vcf"
    params:
        padding = expand("{padding}", padding = config['WES']['PADDING']),
        intervals = expand("{intervals}", intervals = config['WES']['INTERVALS']),
        other = "-G StandardAnnotation -G AS_StandardAnnotation"
    log:
        "logs/gatk_genotype_gvcf/{sample}.log"
    benchmark:
        report("benchmarks/gatk_genotype_gvcf/{sample}.gatkgenotypegvcf", caption = "../report/benchmarking.rst", category = "Benchmarking")
    conda:
        "../envs/gatk4.yaml"
    message:
        "Performing joint genotyping"
    shell:
        """
        gatk --java-options "-Xmx64g -Xms64g" GenotypeGVCFs \
            -R {input.genome} -V {input.vcf} -O {output} {params.padding} {params.intervals} {params.other}
        """