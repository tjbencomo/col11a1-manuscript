# File: helix_breakers.R
# Author: Tomas Bencomo
# Description: Calculate number of helix breaking mutations
# in the triple helical region of COL11A1. A helix breaking
# mutation is any amino acid change  in the triple-helical domain
# that converts a Proline or Glycine to a different amino acid
# that isn't a Proline or Glycine.

library(readr)
library(dplyr)
library(stringr)
library(ggplot2)

# Position indices derived from Uniprot
get_region <- function(position, transcript) {
  # Correction for 38 AA deletion in Isoform C
  if (str_detect(transcript, "ENST00000353414")) {
    if (position > 260) {
      position = position + 38
    }
  }
  # Correction for 115 AA deletion in Isoform 4
  if (str_detect(transcript, "ENST00000512756")) {
    if (position > 299) {
      position = position + 115
    }
  }
  if (position >= 230 & position <= 419) {
    return("Nonhelical region")
  } else if (position >= 420 & position <= 508) {
    return("Triple-helical region (interrupted)")
  } else if (position >= 509 & position <= 511) {
    return("Short nonhelical segment")
  } else if (position >= 512 & position <= 528) {
    return("Telopeptide")
  } else if (position >= 529 & position <= 1542) {
    return("Triple-helical region")
  } else if (position >= 1543 & position <= 1563) {
    return("Nonhelical region (C-terminal)")
  } else {
    return("None")
  }
}

mutations <- read_delim("data/mutations.maf.gz", 
                        "\t", escape_double = FALSE, trim_ws = TRUE, 
                        col_types = c(
                          t_ref_count = "c", t_alt_count = "c", t_depth = "c", 
                          n_depth = "c", n_alt_count = "c", n_ref_count = "c"
                        )) %>%
  mutate(patient = factor(patient)) %>%
  filter(Hugo_Symbol == "COL11A1") 
  
mutations <- mutations %>%
  mutate(ref_aa = str_extract(HGVSp_Short, "[A-Z]"),
         position = str_extract(HGVSp_Short, "[0-9]+"),
         position = as.numeric(position),
         mutant_aa = str_extract(str_sub(HGVSp_Short, 4), "[A-Z]+|\\*"),
         mutant_aa = case_when(
           str_detect(HGVSp_Short, "_splice") ~ "splice",
           str_detect(HGVSp_Short, "=") ~ ref_aa,
           TRUE ~ mutant_aa
         )) %>%
  filter(!is.na(HGVSp_Short))
mutations$region <- apply(mutations[,c('position','Transcript_ID')], 1, 
                          function(x) get_region(as.numeric(x['position']), 
                                                 x['Transcript_ID']))
mutations <- mutations %>%
  mutate(
    pg_mutation = case_when(
      ref_aa %in% c("P", "G") & !(mutant_aa %in%c("P", "G")) ~ 1,
      TRUE  ~ 0
    ),
    helix_breaking = case_when(
      pg_mutation == 1 & region == "Triple-helical region" ~ 1,
      TRUE ~ 0
    )
  )

# mutations <- mutations %>%
#   mutate(helix_breaker = case_when(
#     ref_aa %in% c("P", "G") & !(mutant_aa %in%c("P", "G")) ~ 1,
#     TRUE ~ 0
#   ))

# Plot where mutations fall in the amino acid sequence
p <- mutations %>%
  mutate(pg_mutation = case_when(
    pg_mutation == 0 ~ "Other AA Mutations",
    pg_mutation == 1 ~ "Proline/Glycine Mutation"
  )) %>%
  mutate(pg_mutation = factor(pg_mutation, levels=c("Other AA Mutations", "Proline/Glycine Mutation"))) %>%
  ggplot(aes(position)) +
  geom_histogram(aes(fill = factor(pg_mutation)), bins = 40) +
  facet_grid(rows = vars(pg_mutation)) +
  geom_vline(xintercept = c(529-38, 1542-38), color = "red", linetype = "dashed") + #Note adjustment for Isoform C
  labs(x = "Amino Acid Position in COL11A1", y = "Number Of Mutations") +
  scale_x_continuous(breaks = c(0, 500, 1000, 1500, 1800)) +
  theme(legend.position = "none") +
  annotate("rect", xmin=529-38, xmax=1542-38, ymin=0, ymax=Inf, alpha=0.2, fill="red")
print(p)

# Calculate number of SCCs with at least one helix breaking mutation
# Helix breaking mutation is a P/G replacement in the triple-helical region
prop.samples <- mutations %>%
  group_by(patient, .drop = FALSE) %>%
  summarise(n = sum(helix_breaking)) %>%
  mutate(mutated = case_when(
    n > 0 ~ 1,
    n == 0 ~ 0
  )) %>%
  summarise(prop.broken.helix = sum(mutated) / 100)
print(prop.samples)

print(paste("Total Number of Proline/Glycine Mutations:", sum(mutations$pg_mutation)))
print(paste("Proline/Glycine Mutations in Triple Helix Region:",
            sum(mutations[mutations$region == "Triple-helical region", ]$helix_breaking)))
print(paste("Proline/Glycine Mutations OUTSIDE of Triple Helix Region:",
      sum(mutations[mutations$region != "Triple-helical region", ]$pg_mutation)))
  
