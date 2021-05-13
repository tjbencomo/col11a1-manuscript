# Description: Calculate how many tumors have at least 1 P/G mutation
# in the triple helix domain of either COL2A1, COL11A2, COL5A1, COL5A2,
# or COL11A1. This is to address a reviewer comment about:
# "provide the information or discussion on whether they found any mutations 
# in these collagen alpha genes with similar frequency. Are there any collagen 
# gene mutations that show a strong positive correlation with COL11A1 mutations?"
#
# The transcripts used to annotate mutations are all OK except for COL11A2
# To match transcripts with most expressed and canonical transcripts, I reannotated
# variants using an override for COL11A1 with maf2maf.
#
# NOTE: COL5A2 does not contain a triple helix domain according to UniProt

library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(jaccard)
library(UpSetR)
source("src/helper.R")

fp <- file.path("data", "collagen_annotations.maf")
pfp <- file.path("data", "patient_info.csv")
genes <- c("COL2A1", "COL11A1", "COL11A2", "COL5A1", "COL5A2",
           "COL6A6", "COL22A1", "COL6A3", "COL12A1", "COL14A1")
df <- read_tsv(fp)
patients <- read_csv(pfp)
patients <- patients %>%
  mutate(patient = case_when(
    cohort %in% c("Lee", "South") ~ str_c(patient, ".tumor"),
    TRUE ~ str_c(patient, "-T")
  ))
df$Tumor_Sample_Barcode <- factor(df$Tumor_Sample_Barcode, levels = patients$patient)
df$Hugo_Symbol <- factor(df$Hugo_Symbol, levels = genes)

coding <- c("Missense_Mutation", "Nonsense_Mutation")
muts <- df %>%
  filter(Variant_Classification %in% coding)
muts <- parse_HGVSp(muts)
muts <- annotate_pg_mutations(muts)
muts <- muts %>%
  filter(!is.na(position))

mutmat <- muts %>%
  rowwise() %>% 
  mutate(in_helix = isInHelix(Hugo_Symbol, position)) %>%
  ungroup() %>%
  mutate(helix_breaker = ifelse(pg_mutation & in_helix, 1, 0)) %>%
  group_by(Tumor_Sample_Barcode, Hugo_Symbol, .drop = F) %>%
  summarize(n_muts = sum(helix_breaker)) %>%
  ungroup() %>%
  mutate(mutated = ifelse(n_muts > 0, 1, 0)) %>%
  select(-n_muts) %>%
  pivot_wider(names_from = Hugo_Symbol, values_from = mutated)

print("Number of tumors with a COL2A1 helix breaking mutation:")
print(sum(mutmat$COL2A1))
print("Number of tumors with a COL11A1 helix breaking mutation:")
print(sum(mutmat$COL11A1))
print("Number of tumors with a COL11A2 helix breaking mutation:")
print(sum(mutmat$COL11A2))
print("Number of tumors with a COL5A1 helix breaking mutation:")
print(sum(mutmat$COL5A1))
print("Number of tumors with a COL5A2 helix breaking mutation:")
print(sum(mutmat$COL5A2))

print("Jaccard COL11A1 and COL2A1:")
print(jaccard(mutmat$COL11A1, mutmat$COL2A1))
print("Jaccard COL11A1 and COL11A2:")
print(jaccard(mutmat$COL11A1, mutmat$COL11A2))
print("Jaccard COL11A1 and COL5A1:")
print(jaccard(mutmat$COL11A1, mutmat$COL5A1))

col11a1 <- mutmat$Tumor_Sample_Barcode[mutmat$COL11A1 == 1]
col11a2 <- mutmat$Tumor_Sample_Barcode[mutmat$COL11A2 == 1]
col2a1 <- mutmat$Tumor_Sample_Barcode[mutmat$COL2A1 == 1]
col5a1 <- mutmat$Tumor_Sample_Barcode[mutmat$COL5A1 == 1]
col5a2 <- mutmat$Tumor_Sample_Barcode[mutmat$COL5A2 == 1]

inputList <- list(COL11A1 = col11a1, COL11A2 = col11a2, COL2A1 = col2a1,
                  COL5A1 = col5a1, COL5A2 = col5a2)
upset(fromList(inputList), order.by = "freq")

upset(fromList(inputList))
