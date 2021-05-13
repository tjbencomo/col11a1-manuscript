# File: figure1c.R
# Author: Tomas Bencomo
# Description: Script generates Figure 1C.

library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(reshape2)


get_region_extended <- function(position) {
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


snps <- read_delim("data/mutations.maf.gz", 
                                "\t", escape_double = FALSE, trim_ws = TRUE, 
                                col_types = c(
                                  t_ref_count = "c", t_alt_count = "c", t_depth = "c", 
                                  n_depth = "c", n_alt_count = "c", n_ref_count = "c"
                                ))
coding_types <- c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site")
snps <- snps %>% 
  filter(Hugo_Symbol == "COL11A1", Variant_Classification %in% coding_types)

snps <- snps[complete.cases(snps$HGVSp_Short), ]
extracted_aas <- snps$HGVSp_Short %>% str_extract_all("[A-Z]", simplify = T)
extracted_pos <- snps$HGVSp_Short %>% str_extract("[0-9]+")
snps$reference_aa <- extracted_aas[, 1]
snps$variant_aa <- extracted_aas[, 2]
snps$aa_position <- as.integer(extracted_pos)

snps$region <- sapply(snps$aa_position, get_region_extended)
snps$pg.mutation <- ifelse(snps$reference_aa %in% c("P", "G") & !(snps$variant_aa %in% c("P", "G")), 1, 0)
snps$other.mutation <- ifelse(snps$reference_aa %in% c("P", "G") & !(snps$variant_aa %in% c("P", "G")), 0, 1)


# Note: The triple-helical region's mutation rate is also
# used for the interrupted triple helical region's rate as they
# are both composed of the same domain motif. The triple helical region's
# length takes into account the 89 AAs from the interrupted triple helical region
# to account for this.
region.lengths <- c(229, 190, 89, 3, 17, 1014 + 89, 21, 1806-1563)
region.names <- c("Pre-Nonhelical region", "Nonhelical region", "Triple-helical region (interrupted)",
                  "Short nonhelical segment", "Telopeptide", "Triple-helical region",
                  "Nonhelical region (C-terminal)", "Fibrillar collagen NC1")


snps$region <- factor(snps$region, levels = region.names)


region.info <- data.frame(region = region.names, lengths = region.lengths)
region.info$region <- as.character(region.names)

counts <- snps %>% 
  count(region, .drop = F) %>% 
  inner_join(region.info) %>%
  mutate(normalized.count = n / lengths)

# Count interrupted triple-helical region as triple helical region 
triple.helical.rate <- counts %>% 
  filter(region == "Triple-helical region") %>% 
  pull(normalized.count)
counts <- counts %>%
  mutate(normalized.count = case_when(
    str_detect(region, "Triple-helical") ~ triple.helical.rate,
    TRUE ~ normalized.count
  ))

melted <- melt(select(counts, region, normalized.count), id.var = "region")

p <- ggplot(melted, aes(variable, factor(region, levels = region.names))) + 
  geom_tile(aes(fill = value), width = .4) +
  coord_flip() +
  scale_fill_distiller(palette = "Reds", direction = 1, name = "Mutations/AA", 
                       limits = c(0, .15), breaks = c(0, .05, .1, .15)) +
  theme_classic() +
  labs(x = "", y = "Region")
print(p)

