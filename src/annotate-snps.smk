# File: annotate-snps.smk
# Author: Tomas Bencomo
# Description: Snakemake workflow to annotate variants from
# Durinck and Pickering callsets with VEP. Variants are first
# lifted over to hg38 coordinate space with liftover.py. MAF
# formatted files are then extracted with convert_to_maf.R. 
# MAFs are then re-annotated with maf2maf using VEP and selecting
# variant effects for the default Ensembl transcript.
# Some variants do not match the reference allele and have been excluded
# from analysis. These can be found in the excluded-variants folder

singularity: "docker://continuumio/miniconda3"

ref_fasta = 'data/refs/Homo_sapiens_assembly38.fasta'
vep_dir = 'data/vep_data'
datasets = ['durinck', 'pickering']

wildcard_constraints:
    dataset="|".join(datasets)

rule targets:
    input:
        expand("data/{dataset}.maf", dataset=datasets),

rule maf2maf:
    input:
        maf="data/{dataset}.unannotated.maf",
        fasta=ref_fasta,
        vep_dir=vep_dir
    output:
        "data/{dataset}.maf"
    params:
        tmp_dir=os.path.join("data/{dataset}-tmp")
    conda:
        "../envs/annotation.yml"
    shell:
        """
        mkdir {params.tmp_dir}
        vep_fp=`which vep`
        vep_path=$(dirname "$vep_fp")
        maf2maf.pl --input-maf {input.maf} --output-maf {output} \
            --ref-fasta {input.fasta} --vep-data {input.vep_dir} \
            --vep-path $vep_path --ncbi-build GRCh38 \
            --tmp-dir  {params.tmp_dir} \
            --filter-vcf 0 \
            --cache-version 99
        """


