rule gatk_MarkDuplicates:
    input:
        "../results/mapped/{sample}_sorted.bam"
    output:
        bam = temp("../results/mapped/{sample}_sorted_mkdups.bam"),
        metrics = "../results/mapped/{sample}_sorted_mkdups_metrics.txt"
    params:
        maxmemory = expand('"-Xmx{maxmemory}"', maxmemory = config['MAXMEMORY']),
        tdir = expand("{tdir}", tdir = config['TEMPDIR'])
    log:
        "logs/gatk_MarkDuplicates/{sample}.log"
    benchmark:
        "benchmarks/gatk_MarkDuplicates/{sample}.tsv"
    conda:
        "../envs/gatk4.yaml"
    message:
        "Locating and tagging duplicate reads in {input}"
    shell:
        "gatk MarkDuplicates --java-options {params.maxmemory} -I {input} -O {output.bam} -M {output.metrics} --TMP_DIR {params.tdir} &> {log}"