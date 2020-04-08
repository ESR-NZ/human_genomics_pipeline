rule gatk4_HaplotypeCaller:
    input:
        bams = "mapped/{sample}_bwa_recal.bam"
    output:
        vcf = "vcf/{sample}.raw.snps.indels.AS.g.vcf"
    params:
        genome = expand("{genome}", genome = config["GENOME"]),
        dbsnp = expand("{dbsnp}", dbsnp = config["dbSNP"]),
        tdir = expand("{tdir}", tdir = config["TEMPDIR"]),
        mode = "GVCF"
    log:
        "logs/gatk_haplocall/{sample}.log"
    benchmark:
        "benchmarks/gatk_haplocall/{sample}.gatkhaplocall"
    conda:
        "../envs/gatk4.yaml"
    message:
        "Calling germline SNPs and indels via local re-assembly of haplotypes"
    shell:
        "gatk HaplotypeCaller -I {input.bams} -O {output.vcf} -R {params.genome} -D {params.dbsnp} --tmp-dir {params.tdir} -ERC {params.mode}"