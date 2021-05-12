# File: manuscript_stats.R
# Author: Tomas Bencomo
# Description: This script reproduces statistics
# quoted in the manuscript

library(readr)
library(dplyr)
library(stringr)

load_mutations <- function(filepath) {
  mutations <- read_delim(
    filepath, "\t", escape_double = FALSE, trim_ws = TRUE, 
    col_types = c(
      t_ref_count = "c", t_alt_count = "c", t_depth = "c", 
      n_depth = "c", n_alt_count = "c", n_ref_count = "c"
    )
  ) %>%
    mutate(patient = factor(patient))
}

# 1. Number of SCCs with COL11A1 mutation (Nonsense, Splicing, Missense)
statistic1 <- function(mutations) {
  mutation_types <- c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site")
  stat <- mutations %>%
    filter(
      Hugo_Symbol == "COL11A1", 
      Variant_Classification %in% mutation_types
    ) %>%
    count(patient, .drop = FALSE) %>%
    mutate(mutated = case_when(
      n > 0 ~ 1,
      n == 0 ~ 0
    )) %>%
    summarize(mutated_sccs = sum(mutated)) %>%
    pull()
  output <- "SCCs with at least one Missense, Nonsense, or Splice Site mutation in COL11A1:"
  print(paste(output, stat, "/", length(unique(mutations$patient))))
}

# 2. Percent SCC mutated for top 9 gene hits - see figure1A.R

# 3. Number of SCCs with >= 2 missense mutations
statistic3 <- function(mutations) {
  stat <- mutations %>%
    filter(
      Hugo_Symbol == "COL11A1", 
      Variant_Classification == "Missense_Mutation"
    ) %>%
    count(patient, .drop = FALSE) %>%
    mutate(multiple_mutations = case_when(
      n > 1 ~ 1,
      n <= 1 ~ 0
    )) %>%
    summarize(sccs_multiple_missense = sum(multiple_mutations)) %>%
    pull()
  output <- "SCCs with at least two Missense mutations in COL11A1:"
  print(paste(output, stat))
}

# 4. Number of SCCs with P/G mutation in COL11A1
statistic4 <- function(mutations) {
  stat <- mutations %>%
    filter(Hugo_Symbol == "COL11A1") %>%
    mutate(
      ref_aa = str_extract(HGVSp_Short, "[A-Z]"),
      position = str_extract(HGVSp_Short, "[0-9]+"),
      position = as.numeric(position),
      mutant_aa = str_extract(str_sub(HGVSp_Short, 4), "[A-Z]+|\\*"),
      mutant_aa = case_when(
        str_detect(HGVSp_Short, "_splice") ~ "splice",
        str_detect(HGVSp_Short, "=") ~ ref_aa,
        TRUE ~ mutant_aa
    )) %>%
    mutate(
      pg_mutation = case_when(
        ref_aa %in% c("P", "G") & !(mutant_aa %in%c("P", "G")) ~ 1,
        TRUE  ~ 0
      )
    ) %>%
    filter(pg_mutation == 1) %>%
    count(patient, .drop = FALSE) %>%
    mutate(pg_mutation = case_when(
      n >= 1 ~ 1,
      n == 0 ~ 0
    )) %>%
    summarize(sccs_pg_mutation = sum(pg_mutation)) %>%
    pull()
  output <- "SCCs with at least one Proline or Glycine mutation in COL11A1:"
  print(paste(output, stat))
    
}

# Number of total coding COL11A1 mutations seen cross all samples
statistic5 <- function(mutations) {
  coding_types <- c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site")
  stat <- mutations %>%
    filter(Hugo_Symbol == "COL11A1", Variant_Classification %in% coding_types) %>%
    count() %>%
    pull(n)
  output <- "Total number of COL11A1 mutations: "
  print(paste(output, stat))
}

# Number of COL11A1 mutations that are P/G substitutions (should we include nonsense?)
statistic6 <- function(mutations) {
  coding_types <- c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site")
  stat <- mutations %>%
    filter(Hugo_Symbol == "COL11A1", Variant_Classification %in% coding_types) %>%
    mutate(
      ref_aa = str_extract(HGVSp_Short, "[A-Z]"),
      position = str_extract(HGVSp_Short, "[0-9]+"),
      position = as.numeric(position),
      mutant_aa = str_extract(str_sub(HGVSp_Short, 4), "[A-Z]+|\\*"),
      mutant_aa = case_when(
        str_detect(HGVSp_Short, "_splice") ~ "splice",
        str_detect(HGVSp_Short, "=") ~ ref_aa,
        TRUE ~ mutant_aa
      )) %>%
    mutate(
      pg_mutation = case_when(
        ref_aa %in% c("P", "G") & !(mutant_aa %in%c("P", "G")) ~ 1,
        TRUE  ~ 0
      )) %>%
    summarize(n_pg_muts = sum(pg_mutation)) %>%
    pull(n_pg_muts)
  total_col11a1_coding_muts <- mutations %>%
    filter(Hugo_Symbol == "COL11A1", Variant_Classification %in% coding_types) %>%
    count() %>%
    pull(n)
  output <- "Number of COL11A1 coding mutations that are P/G substitutions: "
  print(paste(output, stat, "/", total_col11a1_coding_muts))
}

# 5. Mutation rates of COL11A1 in TCGA cancers - see Figure 1D's code

# 6. Survival analysis 264 gene signature - see manuscript_survival_analysis.R


data_dir <- "data"
filename <- "mutations.maf.gz"
mutations <- load_mutations(file.path(data_dir, filename))

statistic1(mutations)
statistic3(mutations)
statistic4(mutations)
statistic5(mutations)
statistic6(mutations)
