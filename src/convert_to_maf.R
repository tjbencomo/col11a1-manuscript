# File: convert_to_maf.R 
# Author: Tomas Bencomo
# Description: Convert Cho/Durinck and Pickering mutation callsets that
# have been lifted over to hg38 coordinates into MAFs for re-annotation
# using VEP and maf2maf. This will match the annotation features with
# the South and Lee callsets for easy analysis. See the
# snakemake workflow for reannotation details. Check out the maf2maf
# section of the vcf2maf docs for info on the selected MAF columns.
# Variants without hg38 coordinates are excluded because NA in the position
# field messes up maf2maf.



library(readr)
library(dplyr)
library(stringr)

get_durinck_maf <- function() {
  durinck <- read_csv("data/durinck_mutations_hg38.csv", col_types = c(Chromosome = "c"))
  durinck.maf <- durinck %>% 
    select(
      Chromosome, `Genomic position`, `Genomic reference base`, 
      `Genomic mutated base`, `HGNC symbol`, 
      `Tumor reference base count`, 
      `Tumor mutant base count`, `Normal reference base count`,
      `Normal mutant base count`, patient) %>%
    rename(
      gene = `HGNC symbol`,
      Start_Position = `Genomic position`,
      Reference_Allele = `Genomic reference base`,
      Tumor_Seq_Allele2 = `Genomic mutated base`,
      t_ref_count = `Tumor reference base count`,
      t_alt_count = `Tumor mutant base count`,
      n_ref_count = `Normal reference base count`,
      n_alt_count = `Normal mutant base count`
    ) %>%
    mutate(
      Chromosome = str_c('chr', Chromosome, sep=""),
      t_depth = t_ref_count + t_alt_count,
      n_depth = n_ref_count + n_alt_count,
      Tumor_Sample_Barcode = str_c(patient, "-T", sep=""),
      Matched_Norm_Sample_Barcode = str_c(patient, "-N", sep="")
    ) %>%
    arrange(
      Chromosome,
      Start_Position
    )
}

get_pickering_maf <- function() {
  pickering <- read_csv("data/pickering_mutations_hg38.csv")
  pickering.maf <- pickering %>%
    select(
      Chromosome,
      Start_position,
      End_position,
      Reference_Allele,
      Tumor_Seq_Allele1,
      Tumor_Seq_Allele2,
      Tumor_Sample_Barcode,
      Matched_Norm_Sample_Barcode,
      Hugo_Symbol,
      TTotCov,
      TVarCov,
      NTotCov,
      NVarCov
    ) %>%
    rename(
      gene = Hugo_Symbol,
      t_depth = TTotCov,
      t_alt_count = TVarCov,
      n_depth = NTotCov,
      n_alt_count = NVarCov
    ) %>%
    mutate(
      Chromosome = str_c('chr', Chromosome, sep=""),
      t_depth = as.numeric(na_if(t_depth, "--")),
      n_depth = as.numeric(na_if(n_depth, "--")),
      t_alt_count = as.numeric(na_if(t_alt_count, "--")),
      n_alt_count = as.numeric(na_if(n_alt_count, "--")),
      t_ref_count = t_depth - t_alt_count,
      n_ref_count = n_depth - n_alt_count
    ) %>%
    arrange(
      Chromosome,
      Start_position
    )
}

durinck <- get_durinck_maf()
durinck.missing <- durinck %>%
  filter(is.na(Start_Position))
durinck <- durinck %>%
  tidyr::drop_na(Start_Position)
write_delim(durinck, "data/durinck.unannotated.maf", delim = "\t")
write_delim(durinck.missing, "data/durinck.liftover.problems.maf", delim = "\t")

pickering <- get_pickering_maf()
pickering.missing <- pickering %>%
  filter(is.na(Start_position) | is.na(End_position))
pickering <- pickering %>%
  tidyr::drop_na(Start_position, End_position)
write_delim(pickering, "data/pickering.unannotated.maf", delim = "\t")
write_delim(pickering.missing, "data/pickering.liftover.problems.maf", delim = "\t")

print(str_c(dim(durinck.missing)[1], "variants excluded from Durinck callset because no hg38 coordinates", sep=" "))
print(str_c(dim(pickering.missing)[1], "variants excluded from Pickering callset because no hg38 coordinates", sep = " "))


