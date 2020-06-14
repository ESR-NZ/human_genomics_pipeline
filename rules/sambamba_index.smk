rule sambamba_index:
    input:
        "mapped/{sample}_bwamem_sorted_mkdups.bam"
    output:
        temp("mapped/{sample}_bwamem_sorted_mkdups.bam.bai")
    log:
        "logs/sambamba_index/{sample}.log"
    benchmark:
        "benchmarks/sambamba_index/{sample}.sambamba"
    conda:
        "../envs/sambamba.yaml"
    threads: 8
    message:
        "Building index files for BAM files"
    shell:
        "sambamba index -p {input} -t {threads}"