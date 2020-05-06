rule gatk4_HaplotypeCaller_GVCF:
    input:
        bams = "mapped/{sample}_bwa_recal.bam",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME']),
        dbsnp = expand("{dbsnp}", dbsnp = config['dbSNP'])
    output:
        "vcf/{sample}_haplotype_gvcf.vcf"
    params:
        tdir = expand("{tdir}", tdir = config['TEMPDIR']),
        padding = expand("{padding}", padding = config['WES']['PADDING']),
        intervals = expand("{intervals}", intervals = config['WES']['INTERVALS'])
        other = "-ERC GVCF"
    log:
        "logs/gatk_haplocall/{sample}.log"
    benchmark:
        report("benchmarks/gatk_haplocall/{sample}.gatkhaplocall", caption = "../report/benchmarking.rst", category = "Benchmarking")
    conda:
        "../envs/gatk4.yaml"
    threads: 4
    message:
        "Calling germline SNPs and indels via local re-assembly of haplotypes"
    shell:
        "gatk HaplotypeCaller -I {input.bams} -R {input.refgenome} -D {input.dbsnp} -O {output} --tmp-dir {params.tdir} --native-pair-hmm-threads {threads} {params.padding} {params.intervals} {params.other}"