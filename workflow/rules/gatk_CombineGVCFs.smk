import csv

def get_command(family):
    """Return a string, a portion of the gatk command which defines individuals which should be combined.

    For a particular family, we construct the gatk command by adding -V <individual vcf file> for each individual 
    (defined by individual id column in the pedigree file)
    """
    filename = "../../pedigrees/" + str(family) + "_pedigree.ped"
    
    command = ""
    with open(filename, newline = '') as pedigree:

        pedigree_reader = csv.DictReader(pedigree, fieldnames = ('family', 'individual_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype'), delimiter='\t')
        for individual in pedigree_reader:
            command += "-V ../results/called/" + individual['individual_id'] + "_raw_snps_indels_tmp.g.vcf "

    return command

rule gatk_CombineGVCFs:
    input:
        vcf_dummy = expand("../results/called/{sample}_raw_snps_indels_tmp.g.vcf", sample = SAMPLES), # a dummy vcf to connect this rule to gatk_HaplotypeCaller
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output: 
        vcf = temp("../results/called/{family}_raw_snps_indels_tmp_combined.g.vcf"),
        index = temp("../results/called/{family}_raw_snps_indels_tmp_combined.g.vcf.idx")
    params:
        command = get_command,
        tdir = expand("{tdir}", tdir = config['TEMPDIR']),
        padding = expand("{padding}", padding = config['WES']['PADDING']),
        intervals = expand("{intervals}", intervals = config['WES']['INTERVALS']),
        other = "-G StandardAnnotation -G AS_StandardAnnotation"
    log: 
        "logs/gatk_CombineGVCFs/{family}.log"
    benchmark:
        "benchmarks/gatk_CombineGVCFs/{family}.tsv"
    conda:
        "../envs/gatk4.yaml"
    message:
        "Merging one or more HaplotypeCaller GVCF files into a single GVCF"
    shell:
        "gatk CombineGVCFs -R {input.refgenome} {params.command} -O {output.vcf} --tmp-dir {params.tdir} {params.other} &> {log}"