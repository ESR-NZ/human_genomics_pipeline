human_genomics_pipeline processes paired-end sequencing data (WGS or WES) using `bwa <http://bio-bwa.sourceforge.net/>`_, `sambamba <https://lomereiter.github.io/sambamba/>`_ and `GATK4 <https://gatk.broadinstitute.org/hc/en-us>`_. It is based on the GATK4 best practices for `data pre-processing for variant discovery <https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery>`_

Run parameters:
    * Input data type: {{ snakemake.config["DATA"] }}
    * Reference genome: {{ snakemake.config["REFGENOME"] }}
    * dbSNP database {{ snakemake.config["dbSNP"] }}