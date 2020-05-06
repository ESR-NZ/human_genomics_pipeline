rule sambamba_index:
    input:
        "mapped/{sample}_bwamem_sorted_mkdups.bam"
    output:
        temp("mapped/{sample}_bwamem_sorted_mkdups.bam.bai")
    log:
        "logs/sambamba_index/{sample}.log"
    benchmark:
        report("benchmarks/sambamba_index/{sample}.sambamba", caption = "../report/benchmarking.rst", category = "Benchmarking")
    conda:
        "../envs/sambamba.yaml"
    threads: 4
    message:
        "Building index files for BAM files"
    shell:
        "sambamba index -p {input} -t {threads}"