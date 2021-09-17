rule pbrun_germline:
    input:
        fastq = get_input_fastq,
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output:
        bam = protected("../results/mapped/{sample}_recalibrated.bam"),
        bam_index = protected("../results/mapped/{sample}_recalibrated.bam.bai"),
        vcf = get_output_vcf(config),
        recal = temp("../results/mapped/{sample}_recal.txt")
    resources:
        gpu = config['GPU']
    params:
        readgroup = "--read-group-sm {sample}",
        tdir = expand("{tdir}", tdir = config['TEMPDIR']),
        padding = get_wes_padding_command,
        intervals = get_wes_intervals_command,
        recalibration_resources = get_recal_resources_command,
        other_params = get_params
    log:
        "logs/pbrun_germline/{sample}.log"
    benchmark:
        "benchmarks/pbrun_germline/{sample}.tsv"
    threads: config['THREADS']
    message:
        "Running GPU accelerated germline variant pipeline workflow to generate BAM, vcf and recal output for {input.fastq}"
    shell:
        "/opt/parabricks/3.6/parabricks/pbrun germline --ref {input.refgenome} --in-fq {input.fastq} {params.recalibration_resources} --out-bam {output.bam} --out-variants {output.vcf} --out-recal-file {output.recal} --num-gpus {resources.gpu} {params.readgroup} --tmp-dir {params.tdir} {params.padding} {params.intervals} {params.other_params} --num-cpu-threads {threads} &> {log}"