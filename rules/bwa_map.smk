rule bwa_map:
    input:
        genome = expand("{genome}", genome = config['FILEDIR']['GENOME']),
        R1 = "trim_galore/{sample}_R1_val_1.fq.gz",
        R2 = "trim_galore/{sample}_R2_val_2.fq.gz"
    output: 
        temp("mapped/{sample}_bwamem.bam")
    log:
        "logs/bwamem/{sample}.log"
    benchmark:
        report("benchmarks/bwamem/{sample}.bwamem", caption = "../report/benchmarking.rst", category = "Benchmarking")
    conda:
        "../envs/bwa.yaml"
    threads: 12
    message:
        "Mapping sequences against a reference human genome with BWA-MEM"
    shell: 
        "bwa mem {input.genome} {input.R1} {input.R2} -t {threads} -M | samtools view -@ {threads} -Sbh - > {output}"