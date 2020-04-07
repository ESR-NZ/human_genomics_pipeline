rule sambamba_mkdups:
    input:
        bams = "mapped/{sample}_bwamem_sorted.bam"
    output:
        temp("mapped/{sample}_bwamem_sorted_mkdups.bam")
    params:
        tdir = expand("{tdir}", tdir = config["TEMPDIR"])
    log:
        "logs/sambamba_mkdups/{sample}.log"
    benchmark:
        "benchmarks/sambamba_mkdups/{sample}.sambamba"
    conda:
        "../envs/sambamba.yaml"
    threads: 4
    shell:
        "sambamba markdup -t {threads} {params} --tmpdir={params.tdir} -p {input.bams} {output} --sort-buffer-size=6144 --overflow-list-size=600000 --hash-table-size=600000"