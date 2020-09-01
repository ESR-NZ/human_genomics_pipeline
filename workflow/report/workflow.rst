`human_genomics_pipeline <https://github.com/ESR-NZ/human_genomics_pipeline>`_ processes single samples or cohorts of paired-end sequencing data (WGS or WES) using  `bwa <http://bio-bwa.sourceforge.net/>`_ and `GATK4 <https://gatk.broadinstitute.org/hc/en-us>`_. This workflow is designed to follow the `GATK best practice workflow for germline short variant discovery (SNPs + Indels) <https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels->`_ and is designed to be followed by `vcf_annotation_pipeline <https://github.com/ESR-NZ/vcf_annotation_pipeline>`_.

Run parameters:
    * Input data type: {{ snakemake.config["DATA"] }}
    * Reference genome: {{ snakemake.config["REFGENOME"] }}
    * dbSNP database {{ snakemake.config["dbSNP"] }}
