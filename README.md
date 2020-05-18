# col11a1-manuscript
Code to reproduce analysis and figures for "Increased neoplastic invasion by non-cell autonomous mutant collagen COL11A1" by Lee et. al 2020

## Dependencies
`envs/` has several `.yaml` files that create `conda` environments with the necessary dependencies to run code.
Use `r-env.yaml` to run any of the `.R` scripts. `liftover.yaml` should be used to execute `liftover.py`.
`annotate-snps.smk` is a `snakemake` pipeline. Any python environment with `snakemake` can run the pipeline. It's
recommended to use the `--use-conda` and `--use-singularity` options, which require `conda` and `singularity`. 

## Reproducing Mutation Callset
Follow these steps to reproduce the full mutation callset included with the paper 
(referred to as `mutations.maf.gz` in this repository). 
### 1. Generate mutation callset from Lee lab samples
Mutations were called for samples processed by the Lee lab (53 SCCs in total) using an in-house
variant calling [pipeline](https://github.com/tjbencomo/col11a1-wes-pipeline). The final MAF file
created by the pipeline is referred to as `lee.maf` in this repository. Reads were aligned to hg38.
See the pipeline repository for full details.

### 2. Reannotate Pickering and Cho mutations
Previously published callsets from Pickering (39 SCCs) and Durinck (8 SCCs) were converted to
hg38 coordinates (Pickering was originally in h19 and Durinck hg18) and then reannotated using
Ensembl VEP. 

1. Download the callsets from their respective PubMed links.
2. Run `src/liftover.py` to convert the original callsets from the Pickering and Durinck papers
3. Run `src/convert_to_maf.R` to convert hg38 callsets into MAF format
4. Run `src/annotate-snps.smk` to reannotate the MAFs using `maf2maf`. Note this is a `snakemake` pipeline.
5. Run `src/merge_mafs.R` to combine the Lee, Pickering, and Durinck MAFs into `mutations.maf.gz`.

## Figures
The `src/` folder contains scripts to reproduce the figures. Survival analyses, including Figures 3C 
and Supplemental Figure 5, can be found in `manuscript_survival_analysis.Rmd`
