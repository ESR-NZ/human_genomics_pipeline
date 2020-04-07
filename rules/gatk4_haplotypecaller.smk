rule gatk4_HaplotypeCaller:
    input:
        bams = "mapped/{sample}_bwa_recal.bam"
    output:
        "vcf/{sample}.raw.snps.indels.AS.g.vcf"
    params:
        genome = expand("{genome}", genome = config["GENOME"]),
        dbsnp = expand("{dbsnp}", dbsnp = config["dbSNP"]),
        tdir = expand("{tdir}", tdir = config["TEMPDIR"])
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
        "gatk HaplotypeCaller --reference {params.genome} --emit-ref-confidence GVCF --dbsnp {params.dbsnp} --input {input.bams} --output {output} --tmp-dir {params.tdir}"