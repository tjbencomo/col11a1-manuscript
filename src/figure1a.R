# File: figure1A.R
# Author: Tomas Bencomo
# Description: Script to generate Figure 1A.

library(readr)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)

genes <- c("TP53", "CDKN2A", "COL11A1", "KHDRBS3", "CLASP2",
           "COL4A4", "KNSTRN", "HRAS", "NOTCH1")
mutation_types <- c("Missense_Mutation", "Splice_Site", "Nonsense_Mutation")

patients <- read_csv("data/patient_metadata.csv", col_types = c(patient = "f")) %>%
  mutate(cohort = case_when(
    str_detect(patient, "WD|PD") ~ "South",
    TRUE ~ cohort
  ))
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
         Hugo_Symbol = factor(Hugo_Symbol, levels = rev(genes)),
         Variant_Classification = factor(Variant_Classification))



cohort_order <- c("This study", "South et. al", "Durinck et. al", 
                  "Pickering et. al", "No Mutation")
colors <- c("#97C6AC", "#E3ABB7", "#F4D5A4", "#A4ABC6", "#FFFFFF")

# Compute proportion of SCCs with an "important" mutation
mutated_ratios <- mutations %>%
  count(patient, Hugo_Symbol, .drop = F) %>%
  mutate(mutated = case_when(
    n > 0 ~ 1,
    TRUE ~ 0
  )) %>%
  group_by(Hugo_Symbol) %>%
  summarise(ratio = sum(mutated) / 100) %>%
  mutate(
    percent = str_c(ratio * 100, "%"),
    label = str_c(Hugo_Symbol, " (", percent, ")", sep = "")
    )

label_order <- mutated_ratios$label

# Change fill ordering of cohorts to match sample order
p <- mutations %>%
  count(patient, Hugo_Symbol, .drop = F) %>%
  left_join(patients) %>%
  left_join(mutated_ratios) %>%
  mutate(mutated = case_when(
    cohort == "Lee" & n > 0 ~ "This study",
    cohort == "South" & n > 0 ~ "South et. al",
    cohort == "Pickering" & n > 0 ~ "Pickering et. al",
    cohort == "Durinck" & n > 0 ~ "Durinck et. al",
    TRUE ~ "No Mutation"
  )) %>%
  mutate(id = factor(id),
         mutated = factor(mutated, levels = cohort_order)) %>%
  ggplot(aes(x = id, 
             y = factor(label, levels = label_order), 
             fill = factor(mutated))) +
  geom_tile(colour = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        legend.key = element_rect(colour = 'black'),
        panel.border = element_rect(colour = "black", fill=NA, size = 2),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_fill_manual(values = colors) +
  labs(x = "SCC-Normal Pairs", y = "", fill = "") +
  scale_y_discrete(expand=c(0,0))

print(p)

