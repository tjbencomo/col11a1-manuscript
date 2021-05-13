# Description: Compare nature of COL11A1 mutations between South samples to
# rest of the dataset

library(readr)
library(dplyr)
source("src/helper.R")


fp <- file.path("data", "mutations.maf.gz")
pfp <- file.path("data", "patient_info.csv")
mutations <- read_tsv(fp)
patients <- read_csv(pfp)
patients <- patients %>%
  mutate(patient = case_when(
    cohort %in% c("Lee", "South") ~ str_c(patient, ".tumor"),
    TRUE ~ str_c(patient, "-T")
  ))
mutations$Tumor_Sample_Barcode <- factor(mutations$Tumor_Sample_Barcode, levels = patients$patient)


CODING_TYPES <- c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site")
helix_regions <- c("Triple-helical region", "Triple-helical region (interrupted)")
muts <- mutations %>%
  filter(Hugo_Symbol == "COL11A1", Variant_Classification %in% CODING_TYPES)
muts <- parse_HGVSp(muts)
muts$region <- get_col11a1_regions(muts$position)
muts <- annotate_pg_mutations(muts)

# Chi-square test comparing prevalence of helix-broken tumors in
# South cohort vs other samples
test_number_tumors_with_helix_breaker <- function(muts) {
  mutated <- muts %>%
    mutate(helix_breaker = case_when(
      region %in% helix_regions & pg_mutation ~ 1,
      TRUE ~ 0
    )) %>%
    filter(helix_breaker == 1) %>%
    distinct(Tumor_Sample_Barcode) %>%
    select(Tumor_Sample_Barcode) %>%
    left_join(patients, by = c("Tumor_Sample_Barcode" = "patient")) %>%
    mutate(big_cohort = case_when(
      cohort %in% c("Lee", "Pickering", "Durinck") ~ "Others",
      TRUE ~ "South"
    )) %>%
    count(big_cohort, name = "mutated") %>%
    rename(cohort = big_cohort)
  
  cohort_totals <- data.frame(
    cohort = c("Others", "South"),
    total = c(80, 20)
  )
  
  df <- mutated %>%
    inner_join(cohort_totals) %>%
    mutate(not_mutated = total - mutated) %>%
    select(-cohort, -total) %>%
    as.matrix()
  mat <- t(df)
  print(df)
  print(chisq.test(mat))
  print(fisher.test(mat))
}

# Chi-square test comparing prevalence of Proline/Glycine mutated tumors in
# South cohort vs other samples
test_number_tumors_with_pg_mutation <- function(muts) {
  mutated <- muts %>%
    filter(pg_mutation == 1) %>%
    distinct(Tumor_Sample_Barcode) %>%
    select(Tumor_Sample_Barcode) %>%
    left_join(patients, by = c("Tumor_Sample_Barcode" = "patient")) %>%
    mutate(big_cohort = case_when(
      cohort %in% c("Lee", "Pickering", "Durinck") ~ "Others",
      TRUE ~ "South"
    )) %>%
    count(big_cohort, name = "mutated") %>%
    rename(cohort = big_cohort)
  
  cohort_totals <- data.frame(
    cohort = c("Others", "South"),
    total = c(80, 20)
  )
  
  df <- mutated %>%
    inner_join(cohort_totals) %>%
    mutate(not_mutated = total - mutated) %>%
    select(-cohort, -total) %>%
    as.matrix()
  mat <- t(df)
  print(df)
  print(chisq.test(mat))
  print(fisher.test(mat))
}

# T test comparing mean number of COL11A1 coding mutations in South
# cohort vs mean number of COL11A1 coding mutation in other samples
test_number_col11a1_mutations <- function(muts) {
  df <- muts %>%
    count(Tumor_Sample_Barcode, .drop=FALSE, name = "mutations") %>%
    left_join(patients, by = c("Tumor_Sample_Barcode" = "patient")) %>%
    mutate(collapsed_cohort = case_when(
      cohort %in% c("Lee", "Pickering", "Durinck") ~ "Others",
      TRUE ~ "South"
    )) %>%
    ungroup()
  print(df %>% group_by(collapsed_cohort) %>% summarize(mean_muts = mean(mutations)))
  print(t.test(mutations ~ collapsed_cohort, data=df))
  p <- df %>%
    ggplot(aes(collapsed_cohort, mutations)) +
    geom_boxplot(aes(fill = collapsed_cohort)) +
    theme_bw() +
    labs(x = "Cohort", y = "Number of Mutations In COL11A1") +
    guides(fill = FALSE)
  print(p)
}

# T test comparing mean number of COL11A1 coding mutations in the triple helix region in South
# cohort vs other samples
test_number_col11a1_mutations_in_helix <- function(muts) {
  df <- muts %>%
    filter(region %in% helix_regions) %>%
    count(Tumor_Sample_Barcode, .drop=FALSE, name = "mutations") %>%
    left_join(patients, by = c("Tumor_Sample_Barcode" = "patient")) %>%
    mutate(collapsed_cohort = case_when(
      cohort %in% c("Lee", "Pickering", "Durinck") ~ "Others",
      TRUE ~ "South"
    )) %>%
    ungroup()
  print(df %>% group_by(collapsed_cohort) %>% summarize(mean_muts = mean(mutations)))
  print(t.test(mutations ~ collapsed_cohort, data=df))
  p <- df %>%
    ggplot(aes(collapsed_cohort, mutations)) +
    geom_boxplot(aes(fill = collapsed_cohort)) +
    theme_bw() +
    labs(x = "Cohort", y = "Number of Mutations In Triple Helix") +
    guides(fill = FALSE)
  print(p)
}

test_number_tumors_with_helix_breaker(muts)
test_number_tumors_with_pg_mutation(muts)
test_number_col11a1_mutations(muts)
test_number_col11a1_mutations_in_helix(muts)

