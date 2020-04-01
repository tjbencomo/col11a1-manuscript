# File: figure1b.R
# Author: Tomas Bencomo
# Description: This script creates Figure 1B for the manuscript.


library(readr)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(scales)
library(cowplot)

genes <- c('TP53', 'COL11A1', 'HRAS', 'VWDE', 'KHDRBS3',
           'COL4A4', 'BCLAF1', 'CDKN2A', 'NOTCH1')
mutation_types <- c("Missense_Mutation", "Splice_Site", "Nonsense_Mutation")

patients <- read_csv("data/patient_metadata.csv", col_types = c(patient = "f"))
mutations <- read_delim("data/mutations.maf.gz", 
                        "\t", escape_double = FALSE, trim_ws = TRUE, 
                        col_types = c(
                          t_ref_count = "c", t_alt_count = "c", t_depth = "c", 
                          n_depth = "c", n_alt_count = "c", n_ref_count = "c"
                        )
)  %>%
  filter(Hugo_Symbol %in% genes) %>%
  filter(Variant_Classification %in% mutation_types) %>%
  mutate(patient = factor(patient, levels = patients$patient),
         Hugo_Symbol = factor(Hugo_Symbol, levels = genes),
         Variant_Classification = factor(Variant_Classification))

mutation_order <- c("Splicing", "Nonsense", "Missense")
colors <- c("#407E47", "#5AB1AB", "#80AACD")
mutations %>%
  mutate(Variant_Classification = case_when(
    Variant_Classification == "Splice_Site" ~ "Splicing",
    Variant_Classification == "Nonsense_Mutation" ~ "Nonsense",
    Variant_Classification == "Missense_Mutation" ~ "Missense"
  )) %>%
  mutate(Variant_Classification = factor(Variant_Classification, levels = mutation_order)) %>%
  ggplot(aes(Hugo_Symbol, ..count../sum(..count..))) +
  geom_bar(position = "fill", aes(fill = Variant_Classification)) +
  scale_y_continuous(labels = percent_format(), expand=c(0,0)) +
  labs(x = "", y = "% Mutations") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(colour = 'black', size = 1)) +
  guides(fill=guide_legend(title="Mutation Type")) +
  scale_fill_manual(values = colors)
 
  
