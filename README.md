# col11a1-manuscript
Code for "Mutant collagen COL11A1 enhances invasive cancer" by Lee et. al. 2021

## Data
Mutation data for all 100 SCCs can be downloaded [here](https://drive.google.com/drive/folders/18HBfLd9vzNsC02caiXsDaw7GHDo7VOcx?usp=sharing). 
Save the file to `/data` to ensure scripts run properly without any modification.
The mutations are stored in a gzip compressed MAF formatted file. Any code that refers to `mutations.maf.gz`
is referencing this file. See below for details on how this file was generated.

## Analyses and Figures
The `src/` folder contains scripts to reproduce  analyses and figures. 
Most filenames correspond to their associated figure. Survival analyses, including Figures 3C 
and Supplemental Figure 5, can be found in `manuscript_survival_analysis.Rmd`. 
These scripts can be run using the `r-env` environment (see below about environments).

## Dependencies
`envs/` has several `.yaml` files that create `conda` environments with the necessary dependencies to run code.
If you only wish to recreate figures and analyses, creating the `r-env` environment is sufficent.
If you wish to generate `mutations.maf.gz` from scratch, you will also need to create `col11a1-env`.
### Environments
Before beginning, install Anaconda or Miniconda.

Create the following environments
```
# r-env - used to run scripts that create figures and reproduce analysis
conda env create -f envs/r-env.yaml

# col11a1-env - used to generate full mutation callset from scratch
# only needs to be created if you wish to recreate callset from scratch
conda env create -f envs/liftover.yaml
```

Additionally install `snakemake` in any `conda` environment to run `src/annotate-snps.smk`.
This is only needed if if you wish to recreate the callset from scratch
```
conda install snakemake
```

### Datasets
These datasets are only required to regenerate `mutations.maf.gz` from scratch.

1. Download the hg38 reference FASTA from the Broad's 
[GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle). 
Install `Homo_sapiens_assembly38.fasta` and `Homo_sapiens_assembly38.fasta.fai`.Save the files to `data/refs`

2. Install VEP Cache files by following 
[this](https://uswest.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache) tutorial. 
Download cache data for version 99.
Save the files to `data/vep_data`


## Generating `mutations.maf.gz`
Follow these steps to reproduce the complete mutation callset used in the manuscript. 

### 1. Generate mutation callset from Lee lab samples
Mutations were called for samples processed by the Lee lab (53 SCCs in total) using an in-house
variant calling [pipeline](https://github.com/tjbencomo/col11a1-wes-pipeline). The final MAF file
created by the pipeline is referred to as `lee.maf` in this repository. Reads were aligned to hg38.
See the pipeline repository for details on how to generate `lee.maf`. 

### 2. Reannotate Pickering and Cho mutations
Previously published callsets from [Pickering](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4367811/) (39 SCCs) and [Durinck](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3187561/)
(8 SCCs) were converted to
hg38 coordinates (Pickering was originally in hg19 and Durinck hg18) and then reannotated using
Ensembl VEP. 

1. Download the callsets from their respective PubMed links
2. Inside the `data` directory, create a subdirectory `original-callsets`
3. Place the Pickering and Durinck callsets into `data/original-callsets`. Leave `lee.maf` in `data/`.
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
4. Run `src/annotate-snps.smk` to reannotate the MAFs using `maf2maf`.
```
# Assumes snakemake is installed
# If not run `conda install snakemake` to install
# If singularity is available, it is recommended to include --use-singularity
snakemake -s src/annotate-snps.smk --use-conda
```
5. Run `src/merge_mafs.R` to combine the Lee, Pickering, and Durinck MAFs into `mutations.maf.gz`.
```
conda activate r-env
Rscript src/merge_mafs.R
```
