library(readr)
library(dplyr)
source("src/helper.R")

################################################
## Functions
################################################

total_coding_mutations <- function(mutations) {
  stat <- mutations %>%
    filter(Variant_Classification %in% CODING_TYPES) %>%
    count() %>%
    pull(n)
  return(stat)
}

tumors_with_col11a1_coding_mutation <- function(mutations) {
  stat <- mutations %>%
    filter(Variant_Classification %in% CODING_TYPES) %>%
    distinct(Tumor_Sample_Barcode) %>%
    count() %>%
    pull(n)
  return(stat)
}

tumors_with_twoplus_missense_mutations <- function(mutations) {
  stat <- mutations %>%
    filter(Variant_Classification =="Missense_Mutation") %>%
    count(Tumor_Sample_Barcode, name = "n_muts") %>%
    mutate(mutated = case_when(
      n_muts > 1 ~ "Mutated",
      n_muts <= 1 ~ "Less than 2 mutations"
    )) %>%
    filter(mutated == "Mutated") %>%
    count() %>%
    pull(n)
  return(stat)
}

tumors_with_twoplus_coding_mutations <- function(mutations) {
  stat <- mutations %>%
    filter(Variant_Classification %in% CODING_TYPES) %>%
    count(Tumor_Sample_Barcode, name = "n_muts") %>%
    mutate(mutated = case_when(
      n_muts > 1 ~ "Mutated",
      n_muts <= 1 ~ "Less than 2 mutations"
    )) %>%
    filter(mutated == "Mutated") %>%
    count() %>%
    pull(n)
    return(stat)
}

total_pg_mutations <- function(mutations) {
  stat <- mutations %>%
    filter(Variant_Classification %in% CODING_TYPES) %>%
    summarize(n_pg_mutations = sum(pg_mutation)) %>%
    pull(n_pg_mutations)
  return(stat)
}

total_pg_mutations_in_helix <- function(mutations) {
  helix_regions <- c("Triple-helical region", "Triple-helical region (interrupted)")
  stat <- mutations %>%
    filter(Variant_Classification %in% CODING_TYPES, region %in% helix_regions) %>%
    summarize(n_pg_mutations = sum(pg_mutation)) %>%
    pull(n_pg_mutations)
  return(stat)
}

total_coding_mutations_in_helix <- function(mutations) {
  helix_regions <- c("Triple-helical region", "Triple-helical region (interrupted)")
  stat <- mutations %>%
    filter(Variant_Classification %in% CODING_TYPES, region %in% helix_regions) %>%
    count() %>%
    pull(n)
  return(stat)
}

tumors_with_pg_mutation <- function(mutations) {
  stat <- mutations %>%
    filter(pg_mutation == 1) %>%
    distinct(Tumor_Sample_Barcode) %>%
    count() %>%
    pull(n)
  return(stat)
}

tumors_with_pg_mutation_in_helix <- function(mutations) {
  helix_regions <- c("Triple-helical region", "Triple-helical region (interrupted)")
  stat <- mutations %>%
    filter(pg_mutation == 1, region %in% helix_regions) %>%
    distinct(Tumor_Sample_Barcode) %>%
    count() %>%
    pull(n)
  return(stat)
}

################################################
## ETL
################################################

data_dir <- "data"
mutations_fp <- file.path(data_dir, "mutations.maf.gz")

CODING_TYPES <- c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site")
mutations <- read_tsv(mutations_fp)
col11a1_muts <- mutations %>%
  filter(Hugo_Symbol == "COL11A1")
col11a1_muts <- parse_HGVSp(col11a1_muts)
col11a1_muts$region <- get_col11a1_regions(col11a1_muts$position)
col11a1_muts <- annotate_pg_mutations(col11a1_muts)

################################################
## Statistics
################################################

print("Number of coding COL11A1 mutations:")
print(total_coding_mutations(col11a1_muts))

print("Number of coding COL11A1 mutations in the triple helix:")
print(total_coding_mutations_in_helix(col11a1_muts))

print("Number of proline or glycine mutations in the triple helix:")
print(total_pg_mutations_in_helix(col11a1_muts))

print("Number of COL11A1 mutations that affect a proline or glycine:")
print(total_pg_mutations(col11a1_muts))

print("Number of tumors with 2 or more COL11A1 missense mutations:")
print(tumors_with_twoplus_missense_mutations(col11a1_muts))

print("Number of tumors with one coding mutations in COL11A1:")
print(tumors_with_col11a1_coding_mutation(col11a1_muts))

print("Number of tumors with a COL11A1 mutation affecting a proline or glycine residue:")
print(tumors_with_pg_mutation(col11a1_muts))

print("Number of tumors with a COL11A1 mutation affecting proline or glycine in the triple helix region:")
print(tumors_with_pg_mutation_in_helix(col11a1_muts))

