# File: figureS1A.R
# Author: Tomas Bencomo
# Description:
# Script generates Supplemental Figure 1A that shows ratio of 
# transition vs transversion mutations
# for each SCC and labels them by cohort.

library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(stringr)
library(cowplot)
library(tidyr)
library(patchwork)

# Load mutations
data_dir <- "data"
mutations_fp <- file.path(data_dir, "mutations.maf.gz")
df <- fread(mutations_fp)

# Load patient data
patient_fp <- file.path(data_dir, "patient_info.csv")
patient_info <- fread(patient_fp)

# Define which mtuations are transitions
transitions <- c("A/G", "G/A", "C/T", "T/C")

# Filter mutations and label as transition or transversion
df <- df %>% select(
  patient, 
  VARIANT_CLASS, 
  Reference_Allele, 
  Tumor_Seq_Allele2
  ) %>%
  inner_join(patient_info) %>%
  filter(VARIANT_CLASS == "SNV") %>%
  mutate(
    patient = factor(patient),
    mutation = str_c(Reference_Allele, "/", Tumor_Seq_Allele2)
  ) %>%
  mutate(type = case_when(
    mutation %in% transitions ~ "Transition",
    TRUE ~ "Transversion"
  )) %>%
  mutate(
    type = factor(type, levels = c("Transversion", "Transition"))
  )


# Count number of transitions and transversions per patient
# and transform to percentages
patient_counts <- df %>%
  count(patient, type, .drop = F) %>%
  pivot_wider(names_from = "type", values_from = "n") %>%
  mutate(Total = Transversion + Transition) %>%
  mutate(
    Transversion = Transversion / Total,
    Transition = Transition / Total
  ) %>%
  arrange(desc(Transition)) %>%
  select(-Total) %>%
  pivot_longer(-patient, names_to = "type", values_to = "prop") %>%
  inner_join(patient_info)

# Use to sort by transversion later
sampleid_order <- rev(unique(patient_counts$Sample.ID))

# Define transversion and transition colors
colors <- c("#3D8249", "#00A1D5")

# Define cohort colors
cohort_order <- c("This study", "South et. al", "Durinck et. al", 
                  "Pickering et. al", "No Mutation")
cohort_colors <- c("#97C6AC", "#E3ABB7", "#F4D5A4", "#A4ABC6", "#FFFFFF")



# Sorted by transition percentage
sorted_plot <- patient_counts %>%
  inner_join(patient_info) %>%
  mutate(
    type = factor(type, levels = c("Transversion", "Transition")),
    Sample.ID = factor(Sample.ID, levels = sampleid_order)
  ) %>%
  ggplot(aes(Sample.ID, prop, fill = type)) +
  geom_bar(position = "fill", stat = "identity") + 
  theme_cowplot() +
  guides(fill=guide_legend(title="")) +
  scale_y_continuous(labels = percent_format(), expand=c(0,0), limits = c(0, 1.01)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme(
    axis.text.x = element_text(vjust = 1, angle = 90, size=8),
    legend.text = element_text(face = "plain", color = "black", size = 12)
    ) +
  scale_fill_manual(values = colors) +
  labs(x = "SCC-Normal Pairs", y = "% Mutations")

# Sort by Transition/Transversion percentages
sorted_cohort_plot <- patient_info %>%
  select(Sample.ID, cohort) %>%
  mutate(
    cohort = factor(cohort, levels = c("Lee", "South", "Durinck", "Pickering")),
    Sample.ID = factor(Sample.ID, levels = sampleid_order)
  ) %>%
  pivot_longer(-Sample.ID, names_to = "variable", values_to = "value") %>%
  ggplot(aes(factor(Sample.ID), variable)) +
  geom_tile(aes(fill = value), color = "white") +
  labs(x = "", y = "") +
  guides(fill=guide_legend(title="")) +
  scale_fill_manual(values = cohort_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    legend.text = element_text(face = "plain", color = "black", size = 12),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "cm"),
    legend.position = "right"
  )

# Final figure sorting by transition transverion ratio
sorted_combo <- sorted_cohort_plot / sorted_plot + 
  plot_layout(ncol = 1, heights = c(1, 8))

print(sorted_combo)


## Plot same values but sort SCCs by ID instead of ratio
# # Sorted by SampleID
# unsorted_plot <- patient_counts %>%
#   inner_join(patient_info) %>%
#   mutate(
#     type = factor(type, levels = c("Transversion", "Transition")),
#     Sample.ID = factor(Sample.ID)
#   ) %>%
#   ggplot(aes(Sample.ID, prop, fill = type)) +
#   geom_bar(position = "fill", stat = "identity") + 
#   theme_cowplot() +
#   guides(fill=guide_legend(title="")) +
#   scale_y_continuous(labels = percent_format(), expand=c(0,0), limits = c(0, 1.01)) +
#   scale_x_discrete(expand = c(0, 0)) +
#   theme(
#     axis.text.x = element_text(vjust = 1, angle = 90, size=8),
#     plot.margin = margin(0, 0, 0, 0, "cm")
#   ) +
#   scale_fill_manual(values = colors) +
#   labs(x = "SCC-Normal Pairs", y = "% Mutations")


# # Sorted by SampleID
# unsorted_cohort_plot <- patient_info %>%
#   select(Sample.ID, cohort) %>%
#   mutate(cohort = factor(cohort, levels = c("Lee", "South", "Durinck", "Pickering"))) %>%
#   pivot_longer(-Sample.ID, names_to = "variable", values_to = "value") %>%
#   ggplot(aes(factor(Sample.ID), variable)) +
#   geom_tile(aes(fill = value), color = "white") +
#   labs(x = "", y = "") +
#   guides(fill=guide_legend(title="")) +
#   scale_fill_manual(values = cohort_colors) +
#   theme_bw() +
#   theme(
#     axis.text.x = element_blank(),
#     axis.text.y = element_blank(),
#     axis.ticks = element_blank(),
#     panel.border = element_blank(),
#     panel.grid.major = element_blank(),
#     legend.text = element_text(face = "plain", color = "black", size = 12),
#     panel.grid.minor = element_blank(),
#     plot.margin = margin(0, 0, 0, 0, "cm"),
#     legend.position = "right"
#   )
# 
# # Final figure sorted by SampleID
# unsorted_combo <- unsorted_cohort_plot / unsorted_plot + 
#   plot_layout(ncol = 1, heights = c(1, 8))
