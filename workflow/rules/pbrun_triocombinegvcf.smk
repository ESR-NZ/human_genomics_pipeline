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
            command += "--in-gvcf ../results/called/" + individual['individual_id'] + "_raw_snps_indels_tmp.g.vcf "

    return command

rule pbrun_triocombinegvcf:
    input:
        vcf_dummy = expand("../results/called/{sample}_raw_snps_indels_tmp.g.vcf", sample = SAMPLES), # a dummy vcf to connect this rule to gatk_HaplotypeCaller
        index = expand("../results/called/{sample}_raw_snps_indels_tmp.g.vcf.idx", sample = SAMPLES),
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output:
        temp("../results/called/{family}_raw_snps_indels_tmp_combined.g.vcf")
    resources:
        gpu = config['GPU']
    params:
        command = get_command,
    log:
        "logs/pbrun_triocombinegvcf/{family}.log"
    benchmark:
        "benchmarks/pbrun_triocombinegvcf/{family}.tsv"
    message:
        "Merging one or more HaplotypeCaller GVCF files into a single GVCF"
    shell:
        "pbrun triocombinegvcf --ref {input.refgenome} {params.command} --out-variants {output} &> {log}"