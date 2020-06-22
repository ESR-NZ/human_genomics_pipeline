rule mosdepth:
    input:
        bams = "mapped/{sample}_bwa_recal.bam",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output:
        "coverage/{sample}.per-base.bed.gz",
        "coverage/{sample}.per-base.bed.gz.csi",
        report("coverage/{sample}.mosdepth.summary.txt", caption = "../report/coverage.rst", category = "Coverage"),
        report("coverage/{sample}.mosdepth.global.dist.txt", caption = "../report/coverage.rst", category = "Coverage")
    log:
        "logs/mosdepth/{sample}.log"
    benchmark:
        "benchmarks/mosdepth/{sample}.mosdepth"
    conda:
        "../envs/mosdepth.yaml"
    threads: 16
    message:
        "Undertaking fast BAM/CRAM depth calculation"
    shell: 
        "mosdepth -f {input.refgenome} -t {threads} coverage/{wildcards.sample} {input.bams}"