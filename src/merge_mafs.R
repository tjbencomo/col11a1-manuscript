# File: merge_mafs.R
# Author: Tomas Bencomo
# Description: Merge variants from the three different cohorts: 
# Durinck, Pickering, and Lee. Note Lee includes the South patients.
# The final, combined mutational callset is mutations.maf.gz

library(readr)
library(dplyr)
library(stringr)
library(annotables)

durinck <- read_delim("data/durinck.maf", "\t", escape_double = FALSE, 
                      trim_ws = TRUE, skip = 1,
                      col_types = c(
                        t_alt_count = "c",
                        t_ref_count = "c",
                        t_depth = "c",
                        n_alt_count = "c",
                        n_ref_count = "c",
                        n_depth = "c"
                      )
                    ) %>%
  mutate(Center = "Durinck")
pickering <- read_delim("data/pickering.maf", "\t", escape_double = FALSE, 
                        trim_ws = TRUE, skip = 1, 
                        col_types = c(
                          t_alt_count = "c",
                          t_ref_count = "c",
                          t_depth = "c",
                          n_alt_count = "c",
                          n_ref_count = "c",
                          n_depth = "c"
                          )
                        ) %>%
  mutate(Center = "Pickering")
lee <- read_delim("data/lee.maf", "\t", escape_double = FALSE, trim_ws = TRUE,
                  col_types = c(
                    t_alt_count = "c",
                    t_ref_count = "c",
                    t_depth = "c",
                    n_alt_count = "c",
                    n_ref_count = "c",
                    n_depth = "c"
                  )
                ) %>%
  mutate(Center = "Lee")

durinck <- durinck %>%
  mutate(patient = str_extract(Tumor_Sample_Barcode, "c[A-Z]+[0-9]"))
pickering <- pickering %>%
  mutate(patient = str_extract(Tumor_Sample_Barcode, "CSCC-[0-9]+"))
variants <- lee %>%
  bind_rows(durinck) %>%
  bind_rows(pickering) %>%
  mutate(
    rowid = 1:dim(.)[1],
    chr = str_extract(Chromosome, "[0-9]")
    ) %>%
  left_join(grch38, by=c("Hugo_Symbol" = "symbol", "chr" = "chr")) %>% 
  distinct(rowid, .keep_all = T) %>%
  mutate(
    Entrez_Gene_Id = entrez
  ) %>%
  select(-rowid, -ensgene, -entrez, -chr, -start, -end, 
         -strand, -biotype, -description)

write_delim(variants, "data/mutations.maf.gz", delim = "\t")
