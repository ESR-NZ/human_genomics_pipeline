rule trim_galore_pe:
    input:
        ["../../fastq/{sample}_1.fastq.gz", "../../fastq/{sample}_2.fastq.gz"]
    output:
        temp("../results/trimmed/{sample}_1_val_1.fq.gz"),
        report("../results/trimmed/{sample}_1.fastq.gz_trimming_report.txt", caption = "../report/trimming_R1.rst", category = "Trimming"),
        temp("../results/trimmed/{sample}_2_val_2.fq.gz"),
        report("../results/trimmed/{sample}_2.fastq.gz_trimming_report.txt", caption = "../report/trimming_R1.rst", category = "Trimming")
    params:
        adapters = expand("{adapters}", adapters = config['TRIMMING']['ADAPTERS']),
        other = "-q 20 --paired"
    log:
        "logs/trim_galore/{sample}.log"
    benchmark:
        "benchmarks/trim_galore_pe/{sample}.tsv"
    conda:
        "../envs/trim_galore.yaml"
    threads: 16
    message:
        "Applying quality and adapter trimming to input fastq files: {input}"
    shell:
        "trim_galore {input} -o ../results/trimmed/ {params.adapters} {params.other} -j {threads} &> {log}"
