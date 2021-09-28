## Description: Create table S5b showing correlation between mutations in COL11A1 and other collagens

library(readr)
library(dplyr)

data_dir <- "data"
helix_fp <- file.path(data_dir, "collagen_helix_breaker_correlations.csv")
heterotrimer_fp <- file.path(data_dir, "collagen_heterotrimer_correlation.csv")

helix_corr <- read_csv(helix_fp)
trimer_corr <- read_csv(heterotrimer_fp)

trimer_genes <- c("COL5A1", "COL2A1", "COL5A2", "COL11A2")
helix_genes <- c("COL4A1", "COL19A1", "COL4A3", "COL15A1", "COL4A4")

final_table <- bind_rows(
  helix_corr %>%
    filter(G2 %in% helix_genes),
  trimer_corr %>%
    filter(G2 %in% trimer_genes)
) %>%
  arrange(match(G2, c(trimer_genes, helix_genes)))
