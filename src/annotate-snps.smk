# File: annotate-snps.smk
# Author: Tomas Bencomo
# Description: Snakemake workflow to annotate variants from
# Durinck and Pickering callsets with VEP. Variants are first
# lifted over to hg38 coordinate space with liftover.py. MAF
# formatted files are then extracted with convert_to_maf.R. 
# MAFs are then re-annotated with maf2maf using VEP and selecting
# variant effects for the default Ensembl transcript except for
# COL11A1, which uses a custom transcript specified in override-ensts.
# Some variants do not match the reference allele and have been excluded
# from analysis. These can be found in the excluded-variants folder

singularity: "docker://continuumio/miniconda3"

ref_fasta = '/home/groups/carilee/refs/hg38/Homo_sapiens_assembly38.fasta'
vep_dir = '/home/groups/carilee/refs/vep_data'
datasets = ['durinck', 'pickering']

rule targets:
    input:
        expand("mafs/{dataset}.maf", dataset=datasets),

rule maf2maf:
    input:
        maf="mafs/{dataset}.unannotated.maf",
        fasta=ref_fasta,
        vep_dir=vep_dir
    output:
        "mafs/{dataset}.maf"
    params:
        tmp_dir=os.path.join("{dataset}-tmp")
    conda:
        "envs/annotation.yml"
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


