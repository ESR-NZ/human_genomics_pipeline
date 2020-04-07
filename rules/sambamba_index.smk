rule sambamba_index:
    input:
        bams = "mapped/{sample}_bwamem_sorted_mkdups.bam"
    output:
        temp("mapped/{sample}_bwamem_sorted_mkdups.bam.bai")
    log:
        "logs/sambamba_index/{sample}.log"
    benchmark:
        "benchmarks/sambamba_index/{sample}.sambamba"
    conda:
        "../envs/sambamba.yaml"
    threads: 4
    message:
	 "Building index files for {input.bams}"
    shell:
        "sambamba index -t {threads} -p {input.bams}"