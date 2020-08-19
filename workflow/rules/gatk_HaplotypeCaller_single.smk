rule gatk_HaplotypeCaller_single:
    input:
        bams = "../results/mapped/{sample}_recalibrated.bam",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME']),
        dbsnp = expand("{dbsnp}", dbsnp = config['dbSNP'])
    output:
        protected("../results/called/{sample}_raw_snps_indels.vcf")
    params:
        java = "-Xmx30g",
        tdir = expand("{tdir}", tdir = config['TEMPDIR']),
        padding = expand("{padding}", padding = config['WES']['PADDING']),
        intervals = expand("{intervals}", intervals = config['WES']['INTERVALS'])
    log:
        "logs/gatk_HaplotypeCaller_single/{sample}.log"
    benchmark:
        "benchmarks/gatk_HaplotypeCaller_single/{sample}.tsv"
    conda:
        "../envs/gatk4.yaml"
    threads: 16
    message:
        "Calling germline SNPs and indels via local re-assembly of haplotypes for {input.bams}"
    shell:
        "gatk HaplotypeCaller --java-options {params.java} -I {input.bams} -R {input.refgenome} -D {input.dbsnp} -O {output} --tmp-dir {params.tdir} {params.padding} {params.intervals} --native-pair-hmm-threads {threads} &> {log}"