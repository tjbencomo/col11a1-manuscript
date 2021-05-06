# Description: Create table with most frequent COL11A1 mutations in the triple helix

library(readr)
library(dplyr)
library(stringr)


get_region_element <- function(position) {
  if (position < 230) {
    return("Pre-Nonhelical region")
  } else if (position >= 230 && position <= 419) {
    return("Nonhelical region")
  } else if (position >= 420 && position <= 508) {
    return("Triple-helical region (interrupted)")
    # return("Triple-helical region")
  } else if (position >= 509 && position <= 511) {
    return("Short nonhelical segment")
  } else if (position >= 512 && position <= 528) {
    return("Telopeptide")
  } else if (position >= 529 && position <= 1542) {
    return("Triple-helical region")
  } else if (position >= 1543 && position <= 1563) {
    return("Nonhelical region (C-terminal)")
  } else {
    return("Fibrillar collagen NC1")
  }
}
get_region <- Vectorize(get_region_element)



fp <- file.path("data", "mutations.maf.gz")
muts <- read_tsv(fp)

triple_helix <- c("Triple-helical region", "Triple-helical region (interrupted)")
coding <- c("Missense_Mutation", "Nonsense_Mutation")

df <- muts %>% 
  filter(Hugo_Symbol == "COL11A1", Variant_Classification %in% coding) %>%
  mutate(
    refAA = str_extract(HGVSp_Short, "[A-Z]+"),
    position = as.numeric(str_extract(HGVSp_Short, "[0-9]+")),
    mutAA = str_extract(stringi::stri_reverse(HGVSp_Short), "[A-Z]+|\\*")
  ) %>%
  mutate(
    region = get_region(position)
  ) %>%
  mutate(helix_breaker = case_when(
    region %in% triple_helix & refAA %in% c("P", "G") & !(mutAA %in% c("P", "G")) ~ 1,
    TRUE ~ 0
  ))

tbl <- df %>%
  filter(region %in% triple_helix) %>%
  count(refAA, mutAA, name = "n_mutations") %>%
  arrange(desc(n_mutations)) %>%
  rename(ReferenceAA = "refAA", MutantAA = "mutAA")

write_csv(tbl, "data/top_col11a1_mutations.csv")


