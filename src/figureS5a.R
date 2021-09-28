## Description: Create figure S5A showing mutation correlation with COL11A1

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

col11a1 <- mutmat$patient[mutmat$COL11A1 == 1]
col11a2 <- mutmat$patient[mutmat$COL11A2 == 1]
col2a1 <- mutmat$patient[mutmat$COL2A1 == 1]
col5a1 <- mutmat$patient[mutmat$COL5A1 == 1]
col5a2 <- mutmat$patient[mutmat$COL5A2 == 1]


inputList <- list(COL11A1 = col11a1, COL11A2 = col11a2, COL2A1 = col2a1,
                  COL5A1 = col5a1, COL5A2 = col5a2)


p <- upset(fromList(inputList), order.by = "freq", nsets = length(inputList), nintersects = NA)
print(p)

