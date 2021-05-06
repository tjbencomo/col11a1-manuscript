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

isInHelix <- function(gene, position) {
  if (gene == "COL2A1" & position>= 201 & position <= 1214) {
    return(TRUE)
  } else if (gene == "COL11A2" & position >= 487 & position <= 1500) {
    return(TRUE)
  } else if (gene == "COL5A1" & position >= 559 & position <= 1570) {
    return(TRUE)
  } else if (gene == "COL11A1" & ((position >= 420 & position <= 508) | (position >= 529 & position <= 1542))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

fp <- file.path("data", "collagen_annotations.maf")
pfp <- file.path("data", "patient_info.csv")
df <- read_tsv(fp)
patients <- read_csv(pfp)
patients <- patients %>%
  mutate(patient = case_when(
    cohort %in% c("Lee", "South") ~ str_c(patient, ".tumor"),
    TRUE ~ str_c(patient, "-T")
  ))
df$Tumor_Sample_Barcode <- factor(df$Tumor_Sample_Barcode, levels = patients$patient)

coding <- c("Missense_Mutation", "Nonsense_Mutation")
muts <- df %>%
  filter(Variant_Classification %in% coding) %>%
  mutate(
    refAA = str_extract(HGVSp_Short, "[A-Z]+"),
    position = as.numeric(str_extract(HGVSp_Short, "[0-9]+")),
    mutAA = str_extract(stringi::stri_reverse(HGVSp_Short), "[A-Z]+|\\*")
  ) %>%
  filter(!is.na(position))

# Remove COL5A2 because does not have a triple helix region to be mutated
genes <- c("COL2A1", "COL11A1", "COL11A2", "COL5A1")
res <- muts %>% 
  rowwise() %>% 
  mutate(in_helix = isInHelix(Hugo_Symbol, position))  %>%
  filter(in_helix == 1) %>%
  count(Tumor_Sample_Barcode, Hugo_Symbol, refAA, mutAA) %>%
  filter(refAA %in% c("P", "G") & !(mutAA %in% c("P", "G"))) %>%
  pivot_wider(names_from = "Hugo_Symbol", values_from = "n") %>%
  group_by(Tumor_Sample_Barcode) %>%
  mutate(
    across(everything(), ~replace_na(.x, 0))
  ) %>%
  summarize_at(genes, sum) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(n_mutations = sum(COL2A1, COL11A1, COL11A2, COL5A1)) %>%
  mutate(mutated = ifelse(n_mutations > 0, 1, 0))

print("Number of patients with 1 or more helix breaking mutations in at least one collagen")
print(dim(res)[1])

res2 <- muts %>% 
  rowwise() %>% 
  mutate(in_helix = isInHelix(Hugo_Symbol, position))  %>%
  mutate(helix_breaker = case_when(
    in_helix & refAA %in% c("P", "G") & !(mutAA %in% c("P", "G")) ~ 1,
    TRUE ~ 0
  )) %>%
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
  summarize(num_helix_breakers = sum(helix_breaker)) %>%
  ungroup() %>%
  mutate(broken_helix = ifelse(num_helix_breakers > 0, 1, 0)) %>%
  select(-num_helix_breakers) %>%
  pivot_wider(names_from = "Hugo_Symbol", values_from = "broken_helix") %>%
  group_by(Tumor_Sample_Barcode) %>%
  mutate(
    across(everything(), ~replace_na(.x, 0))
  ) %>%
  # summarize_at(genes, sum) #%>%
  ungroup() %>%
  rowwise() %>%
  mutate(n_mutations = sum(COL2A1, COL11A1, COL11A2, COL5A1)) %>%
  mutate(mutated = ifelse(n_mutations > 0, 1, 0)) -> x

print("Number of patients with 1 or more helix breaking mutations in at least one collagen")
print(sum(res2$mutated))


