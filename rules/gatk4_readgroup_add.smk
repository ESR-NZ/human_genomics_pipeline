rule gatk4_readgroup_add:
    input:
        "mapped/{sample}_bwamem_sorted_mkdups.bam"
    output:
        temp("mapped/{sample}_sorted_mkdups_rgreplaced.bam")
    params:
        tdir = expand("{tdir}", tdir = config["TEMPDIR"]),
        other = "-ID 4 -LB lib1 -PL illumina -PU unit1 -SM 20"
    log:
        "logs/gatk_readgroup/{sample}.log"
    benchmark:
        "benchmarks/gatk_readgroup/{sample}.readgroup"
    conda:
        "../envs/gatk4.yaml"
    message:
        "Assigning all reads to a single new read-group"
    shell:
        """
        gatk AddOrReplaceReadGroups \
        -I {input} \
        -O {output} \
        --TMP_DIR {params.tdir} \
        {params.other}
        """