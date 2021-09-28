# Description: Calculate how many tumors have at least 1 coding mutation
# in list of collagens specified by reviewers/Paul. 
# This is to address a reviewer comment about:
# "provide the information or discussion on whether they found any mutations 
# in these collagen alpha genes with similar frequency. Are there any collagen 
# gene mutations that show a strong positive correlation with COL11A1 mutations?"
#
# The transcripts used to annotate mutations are all OK except for COL11A2/COL6A3/COL14A1
# To match transcripts with most expressed and canonical transcripts, I reannotated
# variants using an override with maf2maf. The correct transcript was determined
# via my skin specific transcript override list
#

library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(jaccard)
library(UpSetR)
source("src/helper.R")

fp <- file.path("data", "collagen_annotations.maf.gz")
pfp <- file.path("data", "patient_info.csv")
# genes <- c("COL2A1", "COL11A1", "COL11A2", "COL5A1", "COL5A2",
#            "COL6A6", "COL22A1", "COL6A3", "COL12A1", "COL14A1")
genes <- c("COL2A1", "COL11A1", "COL11A2", "COL5A1", "COL5A2")

df <- read_tsv(fp)
patients <- read_csv(pfp)
patients <- patients %>%
  mutate(patient = case_when(
    cohort %in% c("Lee", "South") ~ str_c(patient, ".tumor"),
    TRUE ~ str_c(patient, "-T")
  ))
df$patient <- factor(df$Tumor_Sample_Barcode, levels = patients$patient)
df$Hugo_Symbol <- factor(df$Hugo_Symbol, levels = genes)

coding <- c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site")
muts <- df %>%
  filter(Variant_Classification %in% coding)
mutmat <- muts %>%
  count(patient, Hugo_Symbol, .drop = F) %>%
  mutate(isMutated = ifelse(n > 0, 1, 0)) %>%
  select(-n) %>%
  pivot_wider(names_from = Hugo_Symbol, values_from = isMutated)



# muts <- parse_HGVSp(muts)
# muts <- annotate_pg_mutations(muts)
# muts <- muts %>%
#   filter(!is.na(position))

# mutmat <- muts %>%
#   rowwise() %>% 
#   mutate(in_helix = isInHelix(Hugo_Symbol, position)) %>%
#   ungroup() %>%
#   mutate(helix_breaker = ifelse(pg_mutation & in_helix, 1, 0)) %>%
#   group_by(Tumor_Sample_Barcode, Hugo_Symbol, .drop = F) %>%
#   summarize(n_muts = sum(helix_breaker)) %>%
#   ungroup() %>%
#   mutate(mutated = ifelse(n_muts > 0, 1, 0)) %>%
#   select(-n_muts) %>%
#   pivot_wider(names_from = Hugo_Symbol, values_from = mutated)

## Number of coding mutations for each gene
for (g in genes) {
  print(paste("Number of tumors with a", g, "coding mutation:"))
  print(sum(mutmat[[g]]))
}

jc_list <- list()
for (i in 1:length(genes)) {
  print(paste("Jaccard coefficient for COL11A1 and", genes[i], ":"))
  jc <- jaccard(mutmat$COL11A1, mutmat[[genes[i]]])
  print(jc)
  jc_list[[i]] <- jc
}
jcdf <- tibble(G1 = "COL11A1", G2 = genes, jaccard_coef = unlist(jc_list))
print(jcdf)

# print("Number of tumors with a COL2A1 coding mutation:")
# print(sum(mutmat$COL2A1))
# print("Number of tumors with a COL11A1 coding mutation:")
# print(sum(mutmat$COL11A1))
# print("Number of tumors with a COL11A2 coding mutation:")
# print(sum(mutmat$COL11A2))
# print("Number of tumors with a COL5A1 coding mutation:")
# print(sum(mutmat$COL5A1))
# print("Number of tumors with a COL5A2 coding mutation:")
# print(sum(mutmat$COL5A2))



# print("Jaccard COL11A1 and COL2A1:")
# print(jaccard(mutmat$COL11A1, mutmat$COL2A1))
# print("Jaccard COL11A1 and COL11A2:")
# print(jaccard(mutmat$COL11A1, mutmat$COL11A2))
# print("Jaccard COL11A1 and COL5A1:")
# print(jaccard(mutmat$COL11A1, mutmat$COL5A1))

col11a1 <- mutmat$patient[mutmat$COL11A1 == 1]
col11a2 <- mutmat$patient[mutmat$COL11A2 == 1]
col2a1 <- mutmat$patient[mutmat$COL2A1 == 1]
col5a1 <- mutmat$patient[mutmat$COL5A1 == 1]
col5a2 <- mutmat$patient[mutmat$COL5A2 == 1]
# col6a6 <- mutmat$patient[mutmat$COL6A6 == 1]
# col22a1 <- mutmat$patient[mutmat$COL22A1 == 1]
# col6a3 <- mutmat$patient[mutmat$COL6A3 == 1]
# col12a1 <- mutmat$patient[mutmat$COL12A1 == 1]
# col14a1 <- mutmat$patient[mutmat$COL14A1 == 1]

# inputList <- list(COL11A1 = col11a1, COL11A2 = col11a2, COL2A1 = col2a1,
#                   COL5A1 = col5a1, COL5A2 = col5a2, COL6A6 = col6a6,
#                   COL22A1 = col22a1, COL6A3 = col6a3, COL12A1 = col12a1,
#                   COL14A1 = col14a1)

inputList <- list(COL11A1 = col11a1, COL11A2 = col11a2, COL2A1 = col2a1,
                  COL5A1 = col5a1, COL5A2 = col5a2)


upset(fromList(inputList), order.by = "freq", nsets = length(inputList), nintersects = NA)
# upset(fromList(inputList), order.by = "freq", nsets = 6)
# 
# upset(fromList(inputList), nsets = length(inputList))
write_csv(jcdf, "data/collagen_heterotrimer_correlation.csv")
