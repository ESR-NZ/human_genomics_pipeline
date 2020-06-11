rule sambamba_mkdups:
    input:
        "mapped/{sample}_bwamem_sorted.bam"
    output:
        temp("mapped/{sample}_bwamem_sorted_mkdups.bam")
    params:
        tdir = expand("{tdir}", tdir = config['TEMPDIR']),
        other = "--sort-buffer-size=6144 --overflow-list-size=600000 --hash-table-size=600000"
    log:
        "logs/sambamba_mkdups/{sample}.log"
    benchmark:
        "benchmarks/sambamba_mkdups/{sample}.sambamba"
    conda:
        "../envs/sambamba.yaml"
    threads: 4
    message:
        "Finding duplicate reads in BAM files"
    shell:
        "sambamba markdup -p {input} {output} --tmpdir={params.tdir} {params.other} -t {threads}"