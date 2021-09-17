if config['TRIM'] == "No" or config['TRIM'] == "no":
    R1 = "../../fastq/{sample}_1.fastq.gz"
    R2 = "../../fastq/{sample}_2.fastq.gz"
    
if config['TRIM'] == "Yes" or config['TRIM'] == "yes":
    R1 = "../results/trimmed/{sample}_1_val_1.fq.gz"
    R2 = "../results/trimmed/{sample}_2_val_2.fq.gz"

if config['DATA'] == "Single" or config['DATA'] == 'single':
    vcf = "../results/called/{sample}_raw_snps_indels.vcf"
    other_params = ""

if config['DATA'] == "Cohort" or config['DATA'] == 'cohort':
    vcf = "../results/called/{sample}_raw_snps_indels_tmp.g.vcf"
    other_params = "--gvcf"

rule pbrun_germline:
    input:
        R1 = R1,
        R2 = R2,
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output:
        bam = protected("../results/mapped/{sample}_recalibrated.bam"),
        bam_index = protected("../results/mapped/{sample}_recalibrated.bam.bai"),
        vcf = vcf,
        recal = temp("../results/mapped/{sample}_recal.txt")
    resources:
        gpu = config['GPU']
    params:
        readgroup = "--read-group-sm {sample}",
        tdir = expand("{tdir}", tdir = config['TEMPDIR']),
        padding = get_wes_padding_command,
        intervals = get_wes_intervals_command,
        recalibration_resources = get_recal_resources_command,
        other_params = other_params
    log:
        "logs/pbrun_germline/{sample}.log"
    benchmark:
        "benchmarks/pbrun_germline/{sample}.tsv"
    threads: config['THREADS']
    message:
        "Running GPU accelerated germline variant pipeline workflow to generate BAM, vcf and recal output for {input.R1} and {input.R2}"
    shell:
        "pbrun germline --ref {input.refgenome} --in-fq {input.R1} {input.R2} {params.recalibration_resources} --out-bam {output.bam} --out-variants {output.vcf} --out-recal-file {output.recal} --num-gpus {resources.gpu} {params.readgroup} --tmp-dir {params.tdir} {params.padding} {params.intervals} {params.other_params} --num-cpu-threads {threads} &> {log}"