rule sambamba_index_rgadd:
    input:
        bams = "mapped/{sample}_sorted_mkdups_rgreplaced.bam"
    output:
        temp("mapped/{sample}_sorted_mkdups_rgreplaced.bam.bai")
    log:
        "logs/sambamba_index/{sample}_rg.log"
    benchmark:
        "benchmarks/sambamba_index/{sample}_rg.sambamba"
    conda:
        "../envs/sambamba.yaml"
    threads: 4
    message:
	"Building index files for {input.bams}"
    shell:
        "sambamba index -t {threads} -p {input.bams}"