rule pbrun_haplotypecaller_single:
    input:
        bam = "../results/mapped/{sample}_recalibrated.bam",
        index = "../results/mapped/{sample}_recalibrated.bam.bai",
        recal = "../results/mapped/{sample}_recal.txt",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output:
        protected("../results/called/{sample}_raw_snps_indels.vcf")
    resources:
        gpu = config['GPU']
    params:
        tdir = expand("{tdir}", tdir = config['TEMPDIR']),
        padding = expand("{padding}", padding = config['WES']['PADDING']),
        intervals = expand("{intervals}", intervals = config['WES']['INTERVALS'])
    log:
        "logs/pbrun_haplotypecaller/{sample}.log"
    benchmark:
        "benchmarks/pbrun_haplotypecaller/{sample}.tsv"
    message:
        "Calling germline SNPs and indels via local re-assembly of haplotypes for {input.bam}"
    shell:
        "pbrun haplotypecaller --ref {input.refgenome} --in-bam {input.bam} --in-recal-file {input.recal} --out-variants {output} --num-gpus {resources.gpu} --tmp-dir {params.tdir} {params.padding} {params.intervals} &> {log}"