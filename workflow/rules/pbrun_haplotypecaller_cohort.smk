rule pbrun_haplotypecaller_cohort:
    input:
        bam = "../results/mapped/{sample}_recalibrated.bam",
        index = "../results/mapped/{sample}_recalibrated.bam.bai",
        recal = "../results/mapped/{sample}_recal.txt",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output:
        gvcf = temp("../results/called/{sample}_raw_snps_indels_tmp.g.vcf"),
        index = temp("../results/called/{sample}_raw_snps_indels_tmp.g.vcf.idx")
    resources:
        gpu = config['GPU']
    params:
        tdir = expand("{tdir}", tdir = config['TEMPDIR']),
        padding = expand("{padding}", padding = config['WES']['PADDING']),
        intervals = expand("{intervals}", intervals = config['WES']['INTERVALS']),
        other = "--gvcf"
    log:
        "logs/pbrun_haplotypecaller_gvcf/{sample}.log"
    benchmark:
        "benchmarks/pbrun_haplotypecaller_gvcf/{sample}.tsv"
    message:
        "Calling germline SNPs and indels via local re-assembly of haplotypes for {input.bam}"
    shell:
        "pbrun haplotypecaller --ref {input.refgenome} --in-bam {input.bam} --in-recal-file {input.recal} --out-variants {output.gvcf} --num-gpus {resources.gpu} --tmp-dir {params.tdir} {params.padding} {params.intervals} {params.other} &> {log}"