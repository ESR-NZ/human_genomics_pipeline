rule bwa_mem:
    input:
        fastq = get_input_fastq,
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output: 
        temp("../results/mapped/{sample}_sorted.bam")
    params:
        readgroup = "'@RG\\tID:{sample}_rg1\\tLB:lib1\\tPL:bar\\tSM:{sample}\\tPU:{sample}_rg1'",
        sortsam = "--MAX_RECORDS_IN_RAM=5000000 --SORT_ORDER=coordinate -I=/dev/stdin",
        tdir = expand("{tdir}", tdir = config['TEMPDIR'])
    log:
        "logs/bwa_mem/{sample}.log"
    benchmark:
        "benchmarks/bwa_mem/{sample}.tsv"
    conda:
        "../envs/bwa.yaml"
    threads: config['THREADS']
    resources:
        mem_mb = config['MAXMEMORY'],
        partition = config['PARTITION']['CPU']
    message:
        "Mapping sequences against a reference human genome with BWA-MEM for {input.fastq}"
    shell:
        'bwa mem '
        '-R {params.readgroup} '
        '-t {threads} '
        '-K 10000000 {input.refgenome} {input.fastq} | '
        'gatk SortSam '
        '--java-options "-Xmx{resources.mem_mb}m" '
        '{params.sortsam} '
        '-O={output} '
        '--TMP_DIR={params.tdir} '
        '&> {log}'