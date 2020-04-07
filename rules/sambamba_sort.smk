rule sambamba_sort:
    input:
        bams = "mapped/{sample}_bwamem.bam"
    output:
        temp("mapped/{sample}_bwamem_sorted.bam")
    params:
        tdir = expand("{tdir}", tdir = config["TEMPDIR"])
    log:
        "logs/sambamba_sort/{sample}.log"
    benchmark:
        "benchmarks/sambamba_sort/{sample}.sambamba"
    conda:
        "../envs/sambamba.yaml"
    threads: 4
    message:
	    "Sorting BAM files"
    shell:
        "sambamba sort -t {threads} -m 6G --tmpdir={params.tdir} -p -o {output} {input.bams}"