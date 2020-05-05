rule gatk4_HaplotypeCaller_GVCF:
    input:
        bams = "mapped/{sample}_bwa_recal.bam",
        genome = expand("{genome}", genome = config["FILEDIR"]["GENOME"]),
        dbsnp = expand("{dbsnp}", dbsnp = config["FILEDIR"]["dbSNP"])
    output:
        "vcf/{sample}_haplotype_gvcf.vcf"
    params:
        tdir = expand("{tdir}", tdir = config["TEMPDIR"]),
        other = "-ERC GVCF"
    log:
        "logs/gatk_haplocall/{sample}.log"
    benchmark:
        "benchmarks/gatk_haplocall/{sample}.gatkhaplocall"
    conda:
        "../envs/gatk4.yaml"
    threads: 4
    message:
        "Calling germline SNPs and indels via local re-assembly of haplotypes"
    shell:
        "gatk HaplotypeCaller -I {input.bams} -R {input.genome} -D {input.dbsnp} -O {output} --tmp-dir {params.tdir} {params.other} --native-pair-hmm-threads {threads}"