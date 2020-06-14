rule sambamba_index_rgadd:
    input:
        "mapped/{sample}_sorted_mkdups_rgreplaced.bam"
    output:
        temp("mapped/{sample}_sorted_mkdups_rgreplaced.bam.bai")
    log:
        "logs/sambamba_index/{sample}_rg.log"
    benchmark:
        "benchmarks/sambamba_index/{sample}_rg.sambamba"
    conda:
        "../envs/sambamba.yaml"
    threads: 8
    message:
        "Building index files for BAM files"
    shell:
        "sambamba index -p {input} -t {threads}"