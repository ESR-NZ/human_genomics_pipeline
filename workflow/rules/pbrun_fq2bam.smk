rule pbrun_fq2bam:
    input:
        R1 = "../results/trimmed/{sample}_1_val_1.fq.gz",
        R2 = "../results/trimmed/{sample}_2_val_2.fq.gz",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output:
        bam = protected("../results/mapped/{sample}_recalibrated.bam"),
        index = protected("../results/mapped/{sample}_recalibrated.bam.bai"),
        recal = temp("../results/mapped/{sample}_recal.txt")
    resources:
        gpu = 1
    params:
        readgroup = "--read-group-sm {sample}",
        tdir = expand("{tdir}", tdir = config['TEMPDIR']),
        padding = expand("{padding}", padding = config['WES']['PADDING']),
        intervals = expand("{intervals}", intervals = config['WES']['INTERVALS']),
        recalibration_resources = expand("{recalibration_resources}", recalibration_resources = config['RECALIBRATION']['RESOURCES'])
    log:
        "logs/pbrun_fq2bam/{sample}.log"
    benchmark:
        "benchmarks/pbrun_fq2bam/{sample}.tsv"
    message:
        "Generating a BAM output for {input.R1} and {input.R2} using BWA-Mem, gatk MarkDuplicates and gatk BaseRecalibrator"
    shell:
        "pbrun fq2bam --ref {input.refgenome} --in-fq {input.R1} {input.R2} {params.recalibration_resources} --out-bam {output.bam} --out-recal {output.recal} {params.readgroup} --tmp-dir {params.tdir} {params.padding} {params.intervals} > {log}"