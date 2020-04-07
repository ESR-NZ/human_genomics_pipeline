rule bwa_map:
    input:
        R1 = "trim_galore/{sample}_R1_val_1.fq.gz",
        R2 = "trim_galore/{sample}_R2_val_2.fq.gz"
    output: 
        temp("mapped/{sample}_bwamem.bam")
    params:
        genome = expand("{genome}", genome = config["GENOME"])
    log:
        "logs/bwamem/{sample}.log"
    benchmark:
        "benchmarks/bwamem/{sample}.bwamem"
    conda:
        "../envs/bwa.yaml"
    threads: 12
    shell: 
        "bwa mem -M -t {threads} {params.genome} {input.R1} {input.R2} | samtools view -@ {threads} -Sbh - > {output}"