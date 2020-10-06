if config['TRIM'] == "No" or config['TRIM'] == "no":
    R1 = "../../fastq/{sample}_1.fastq.gz"
    R2 = "../../fastq/{sample}_2.fastq.gz"
    
if config['TRIM'] == "Yes" or config['TRIM'] == "yes":
    R1 = "../results/trimmed/{sample}_1_val_1.fq.gz"
    R2 = "../results/trimmed/{sample}_2_val_2.fq.gz"

rule bwa_mem:
    input:
        R1 = R1,
        R2 = R2,
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output: 
        temp("../results/mapped/{sample}_sorted.bam")
    params:
        readgroup = "'@RG\\tID:{sample}_rg1\\tLB:lib1\\tPL:bar\\tSM:{sample}\\tPU:{sample}_rg1'",
        maxmemory = expand('"-Xmx{maxmemory}"', maxmemory = config['MAXMEMORY']),
        sortsam = "--MAX_RECORDS_IN_RAM=5000000 --SORT_ORDER=coordinate -I=/dev/stdin",
        tdir = expand("{tdir}", tdir = config['TEMPDIR'])
    log:
        "logs/bwa_mem/{sample}.log"
    benchmark:
        "benchmarks/bwa_mem/{sample}.tsv"
    conda:
        "../envs/bwa.yaml"
    threads: config['THREADS']
    message:
        "Mapping sequences against a reference human genome with BWA-MEM for {input.R1} {input.R2}"
    shell:
        "bwa mem -R {params.readgroup} -t {threads} -K 10000000 {input.refgenome} {input.R1} {input.R2} | gatk SortSam --java-options {params.maxmemory} {params.sortsam} -O={output} --TMP_DIR={params.tdir} &> {log}"