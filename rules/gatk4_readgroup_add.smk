rule gatk4_readgroup_add:
    input:
        bams = "mapped/{sample}_bwamem_sorted_mkdups.bam"
    output:
        temp("mapped/{sample}_sorted_mkdups_rgreplaced.bam")
    log:
        "logs/gatk_readgroup/{sample}.log"
    benchmark:
        "benchmarks/gatk_readgroup/{sample}.readgroup"
    conda:
        "../envs/gatk4.yaml"
    threads: 4
    shell:
        "gatk AddOrReplaceReadGroups --INPUT {input.bams} --OUTPUT {output} --RGID 4 --RGLB lib1 --RGPL illumina --RGPU unit1 --RGSM 20"