# col11a1-manuscript
Code for "Increased neoplastic invasion by non-cell autonomous mutant collagen COL11A1" by Lee et. al. 2020

## Dependencies
`envs/` has several `.yaml` files that create `conda` environments with the necessary dependencies to run code.
Use `r-env.yaml` to run any of the `.R` scripts. `liftover.yaml` should be used to execute `liftover.py`.
`annotate-snps.smk` is a `snakemake` pipeline. Any python environment with `snakemake` can run the pipeline. It's
recommended to use the `--use-conda` and `--use-singularity` options, which require `conda` and `singularity`. 
### Environments
Before beginning, install Anaconda or Miniconda.


Create the following environments
```
# col11a1-env
conda env create -f envs/liftover.yaml

# r-env
conda env create -f envs/r-env.yaml
```

Additionally install `snakemake` a `conda` environment of your choice to run `src/annotate-snps.smk`
```
conda install snakemake
```

### Datasets
1. Download the hg38 reference FASTA from the Broad's 
[GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle). 
Install `Homo_sapiens_assembly38.fasta` and `Homo_sapiens_assembly38.fasta.fai`.Save the files to `data/refs`

2. Install VEP Cache files by following 
[this](https://uswest.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache) tutorial. 
Download cache data for version 99.
Save the files to `data/vep_data`


## Generating the Full Mutation Callset
Follow these steps to reproduce the full mutation callset used for analysis in the paper 
(referred to as `mutations.maf.gz` in this repository).

### 1. Generate mutation callset from Lee lab samples
Mutations were called for samples processed by the Lee lab (53 SCCs in total) using an in-house
variant calling [pipeline](https://github.com/tjbencomo/col11a1-wes-pipeline). The final MAF file
created by the pipeline is referred to as `lee.maf` in this repository. Reads were aligned to hg38.
See the pipeline repository for full details. 

The fully annotated mutations for the 53 SCCs we processed are also available to download from our manuscript
on the journal's website. Download the supplementary file and rename it `lee.maf` if you'd prefer
to skip the pipeline. 

### 2. Reannotate Pickering and Cho mutations
Previously published callsets from [Pickering](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4367811/) (39 SCCs) and [Durinck](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3187561/)
(8 SCCs) were converted to
hg38 coordinates (Pickering was originally in hg19 and Durinck hg18) and then reannotated using
Ensembl VEP. 

1. Download the callsets from their respective PubMed links
2. Inside the `data` directory, create a subdirectory `original-callsets`
3. Place the Pickering, Durinck, and Lee callsets into `data/original-callsets`
2. Run `src/liftover.py` to convert the original callsets from the Pickering and Durinck papers
```
conda activate col11a1-env
python src/liftover.py
```
3. Run `src/convert_to_maf.R` to convert hg38 callsets into MAF format
```
conda activate r-env
Rscript src/convert_to_maf.R
```
4. Run `src/annotate-snps.smk` to reannotate the MAFs using `maf2maf`. Note this is a `snakemake` pipeline.
```
# Assumes snakemake is installed
# If not run `conda install snakemake` to install
# If singularity is available, it is recommended to use --use-singularity
snakemake -s src/annotate-snps.smk --use-conda
```
5. Run `src/merge_mafs.R` to combine the Lee, Pickering, and Durinck MAFs into `mutations.maf.gz`.
```
conda activate r-env
Rscript src/merge_mafs.R
```

## Figures
The `src/` folder contains scripts to reproduce the figures. 
Most filenames correspond to their associated figure. Survival analyses, including Figures 3C 
and Supplemental Figure 5, can be found in `manuscript_survival_analysis.Rmd`. 
These scripts can be run using the `r-env` environment.
