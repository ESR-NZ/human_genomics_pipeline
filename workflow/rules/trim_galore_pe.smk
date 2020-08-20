rule trim_galore_pe:
    input:
        ["../../fastq/{sample}_1.fastq.gz", "../../fastq/{sample}_2.fastq.gz"]
    output:
        temp("../results/trimmed/{sample}_1_val_1.fq.gz"),
        report("../results/trimmed/{sample}_1.fastq.gz_trimming_report.txt", caption = "../report/trimming_R1.rst", category = "Trimming"),
        temp("../results/trimmed/{sample}_2_val_2.fq.gz"),
        report("../results/trimmed/{sample}_2.fastq.gz_trimming_report.txt", caption = "../report/trimming_R1.rst", category = "Trimming")
    params:
        extra = "--illumina -q 20"
    log:
        "logs/trim_galore/{sample}.log"
    benchmark:
        "benchmarks/trim_galore_pe/{sample}.tsv"
    message:
        "Applying quality and adapter trimming of input fastq files: {input}"
    wrapper:
        "0.64.0/bio/trim_galore/pe"