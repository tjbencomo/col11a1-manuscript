library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(jaccard)
source("src/helper.R")

genes <- c("COL19A1", "COL6A3", "COL4A4", "COL6A6", "COL15A1", 
           "COL3A1", "COL4A1", "COL4A3", "COL6A5", "COL4A2", 
           "COL11A1")
coding <- c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site")
fp <- "data/mutations.maf.gz"
pfp <- file.path("data", "patient_info.csv")
df <- read_tsv(fp)
patients <- read_csv(pfp)
patients <- patients %>%
  mutate(patient = case_when(
    cohort %in% c("Lee", "South") ~ str_c(patient, ".tumor"),
    TRUE ~ str_c(patient, "-T")
  ))
df <- df %>%
  filter(Hugo_Symbol %in% genes, Variant_Classification %in% coding)
df$patient <- factor(df$Tumor_Sample_Barcode, levels = patients$patient)
df$Hugo_Symbol <- factor(df$Hugo_Symbol, levels = genes)



muts <- parse_HGVSp(df)
muts <- annotate_pg_mutations(muts)
muts <- muts %>%
  filter(!is.na(position))

mutmat <- muts %>%
  rowwise() %>%
  mutate(in_helix = isInHelix(Hugo_Symbol, position)) %>%
  ungroup() %>%
  mutate(helix_breaker = ifelse(pg_mutation & in_helix, 1, 0)) %>%
  group_by(patient, Hugo_Symbol, .drop = F) %>%
  summarize(n_muts = sum(helix_breaker)) %>%
  ungroup() %>%
  mutate(mutated = ifelse(n_muts > 0, 1, 0)) %>%
  select(-n_muts) %>%
  pivot_wider(names_from = Hugo_Symbol, values_from = mutated)

jc_list <- list()
for (i in 1:length(genes)) {
  print(paste("Jaccard coefficient for COL11A1 and", genes[i], ":"))
  jc <- jaccard(mutmat$COL11A1, mutmat[[genes[i]]])
  print(jc)
  jc_list[[i]] <- jc
}
jcdf <- tibble(G1 = "COL11A1", G2 = genes, jaccard_coef = unlist(jc_list))
print(jcdf)
