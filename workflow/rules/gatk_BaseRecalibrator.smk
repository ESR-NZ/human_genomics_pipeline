def get_recal_resources_command(resource):
    """Return a string, a portion of the gatk BaseRecalibrator command (used in the gatk_BaseRecalibrator and the
    parabricks_germline rules) which dynamically includes each of the recalibration resources defined by the user
    in the configuration file. For each recalibration resource (element in the list), we construct the command by
    adding either --knownSites (for parabricks) or --known-sites (for gatk4) <recalibration resource file>
    """
    
    recal_command = ""
    
    for resource in config['RECALIBRATION']['RESOURCES']:
        if config['GPU_ACCELERATED'] == "Yes" or config['GPU_ACCELERATED'] == "yes":
            recal_command += "--knownSites " + resource + " "

        if config['GPU_ACCELERATED'] == "No" or config['GPU_ACCELERATED'] == "no":
            recal_command += "--known-sites " + resource + " "

    return recal_command

rule gatk_BaseRecalibrator:
    input:
        bams = "../results/mapped/{sample}_sorted_mkdups.bam",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output:
        report("../results/mapped/{sample}_recalibration_report.grp", caption = "../report/recalibration.rst", category = "Base recalibration")
    params:
        maxmemory = expand('"-Xmx{maxmemory}"', maxmemory = config['MAXMEMORY']),
        tdir = expand("{tdir}", tdir = config['TEMPDIR']),
        padding = expand("{padding}", padding = config['WES']['PADDING']),
        intervals = expand("{intervals}", intervals = config['WES']['INTERVALS']),
        recalibration_resources = get_recal_resources_command
    log:
        "logs/gatk_BaseRecalibrator/{sample}.log"
    benchmark:
        "benchmarks/gatk_BaseRecalibrator/{sample}.tsv"
    conda:
        "../envs/gatk4.yaml"
    message:
        "Generating a recalibration table for {input.bams}"
    shell:
        "gatk BaseRecalibrator --java-options {params.maxmemory} -I {input.bams} -R {input.refgenome} -O {output} --tmp-dir {params.tdir} {params.padding} {params.intervals} {params.recalibration_resources} &> {log}"