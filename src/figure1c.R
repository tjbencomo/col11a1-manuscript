# File: figure1c.R
# Author: Tomas Bencomo
# Description: Script generates Figure 1C. This script creates
# a heatmap for COL11A1 Protein Isoform C. Isoform C differs from Isoform A by a deletion
# from AA 261 - 299. While this deletion affects a nonhelical region,
# the indexing for the other AAs in the protein are shifted leftward by 38 AAs
# after the 260th AA position. get_region_extended() accounts for this shift.

library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(reshape2)


get_region_extended <- function(position) {
  if (position > 260) {
    position = position + 38
  }
  if (position < 230) {
    return("Pre-Nonhelical region")
  } else if (position >= 230 && position <= 419) {
    return("Triple-helical region (interrupted)")
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


snps <- read_delim("data/mutations.maf.gz", 
                                "\t", escape_double = FALSE, trim_ws = TRUE, 
                                col_types = c(
                                  t_ref_count = "c", t_alt_count = "c", t_depth = "c", 
                                  n_depth = "c", n_alt_count = "c", n_ref_count = "c"
                                ))
snps <- snps %>% filter(Hugo_Symbol == "COL11A1", Variant_Classification %in% c("Missense_Mutation",
                                                                                "Nonsense_Mutation"))

snps <- snps[complete.cases(snps$HGVSp_Short), ]
extracted_aas <- snps$HGVSp_Short %>% str_extract_all("[A-Z]", simplify = T)
extracted_pos <- snps$HGVSp_Short %>% str_extract("[0-9]+")
snps$reference_aa <- extracted_aas[, 1]
snps$variant_aa <- extracted_aas[, 2]
snps$aa_position <- as.integer(extracted_pos)

snps$region <- sapply(snps$aa_position, get_region_extended)
snps$pg.mutation <- ifelse(snps$reference_aa %in% c("P", "G") & !(snps$variant_aa %in% c("P", "G")), 1, 0)
snps$other.mutation <- ifelse(snps$reference_aa %in% c("P", "G") & !(snps$variant_aa %in% c("P", "G")), 0, 1)

region.lengths <- c(229, 419 - 230, 511-509, 528-512, 1542-529, 1563-1543, 1767-1564)
region.names <- c("Pre-Nonhelical region", "Triple-helical region (interrupted)",
                  "Short nonhelical segment", "Telopeptide", "Triple-helical region",
                  "Nonhelical region (C-terminal)", "Fibrillar collagen NC1")

snps$region <- factor(snps$region, levels = region.names)


region.info <- data.frame(region = region.names, lengths = region.lengths)
region.info$region <- as.character(region.names)

counts <- snps %>% 
  count(region, .drop = F) %>% 
  inner_join(region.info) %>%
  mutate(normalized.count = n / lengths)

melted <- melt(select(counts, region, normalized.count), id.var = "region")

p <- ggplot(melted, aes(variable, factor(region, levels = region.names))) + 
  geom_tile(aes(fill = value), width = .4) +
  coord_flip() +
  scale_fill_distiller(palette = "Reds", direction = 1, name = "Mutations/AA", 
                       limits = c(0, .15), breaks = c(0, .05, .1, .15)) +
  labs(x = "", y = "Region")
print(p)

