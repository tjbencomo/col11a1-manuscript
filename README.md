# col11a1-manuscript
Code to reproduce figures for "Increased neoplastic invasion by non-cell autonomous mutant collagen COL11A1" by Lee et. al 2020

## Notes
SCCs sequenced by the Lee lab were analyzed using an in-house somatic variant calling [pipeline](https://github.com/tjbencomo/col11a1-wes-pipeline). The pipeline created `lee.maf` by aggregating variants
from each sample and annotating them using VEP and vcf2maf. Reads were aligned to hg38. 

Callsets from Durinck and Pickering were converted to hg38 coordinates. 
Durinck originally used hg18 coordinates. 
Pickering used hg19 coordinates. See `liftover.py` for details. 
The converted Durinck and Pickering callsets were transformed into MAFs via `convert_to_maf.R` and
then reannotated using VEP and maf2maf. See the Snakemake workflow `annotate-snps.smk` for more details. 
`merge_mafs.R` aggregates the three annotated MAFs into a single
combined MAF, `mutations.maf.gz`. 
This file was used for analysis and figure generation.

## Dependencies
`liftover.py` requires the `pyliftover` library. You can build a conda environment with `envs/liftover.yaml`
to run `liftover.py`
