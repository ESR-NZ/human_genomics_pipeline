rule gatk4_readgroup_add:
    input:
        bams = "mapped/{sample}_bwamem_sorted_mkdups.bam"
    output:
        temp("mapped/{sample}_sorted_mkdups_rgreplaced.bam")
    params:
        extra = "-ID 4 -LB lib1 -PL illumina -PU unit1 -SM 20"
    log:
        "logs/gatk_readgroup/{sample}.log"
    benchmark:
        "benchmarks/gatk_readgroup/{sample}.readgroup"
    conda:
        "../envs/gatk4.yaml"
    threads: 4
    message:
	    "Assigning all reads to a single new read-group"
    shell:
        "gatk AddOrReplaceReadGroups -I {input.bams} -O {output} {params.extra}"