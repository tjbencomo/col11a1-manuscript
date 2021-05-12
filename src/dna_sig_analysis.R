# Description: Figure S1A shows the majority of SCCs from South's cohort
# have mostly transversion mutations and very few transition mutations. This
# is in contrast to to the other samples from the other cohorts, which for the
# most part are majority transition (transition mutations arise from UV damage). 
# South's cohort is mostly made up of immunocompromised individuals, and we suspect
# this mutation signature difference is due to the immunocompromised (IC) status of
# these patients. This script investigates this phenomena.

library(readr)
library(data.table)
library(cowplot)
library(readxl)
library(scales)
library(tidyr)
library(stringr)
library(dplyr)
library(ggplot2)
library(forcats)

# Load mutations
data_dir <- "data"
mutations_fp <- file.path(data_dir, "mutations.maf.gz")
south_meta_fp <- file.path(data_dir, "NIHMS57574-supplement-data_tables.xlsx")

df <- fread(mutations_fp)

# Load patient data
patient_fp <- file.path(data_dir, "patient_info.csv")
patient_info <- fread(patient_fp)
patient_info <- patient_info %>% 
  mutate(
    cohort = case_when(
      str_detect(patient, "WD|PD") ~ "South",
      TRUE ~ "Other"
    ),
    oid = as.numeric(str_extract(patient, "[0-9]+")),
    diff_status = case_when(
      str_detect(patient, "WD") ~ "WD",
      str_detect(patient, "PD") ~ "PD",
      TRUE ~ "Unknown"
    )
  )

## Load south metadata
south_met <- read_excel(south_meta_fp, skip = 3)
south_met <- south_met %>%
  mutate(
    cohort = "South",
    oid = as.numeric(str_extract(Patient, "[0-9]+")),
    diff_status = str_extract(Patient, "[A-Z]+")
  ) %>%
  rename(ImmuneStatus = `Immune Status`)

patient_info <- patient_info %>%
  left_join(south_met, by = c("cohort", "oid", "diff_status"))

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

#########################################################
## Cohort Comparison
#########################################################


#Compare transversion level in PD vs other samples
with(
  patient_counts %>%
    filter(diff_status != "WD", type == "Transversion") %>%
    mutate(dstatus = ifelse(diff_status == "PD", "PD", "Other")),
  wilcox.test(prop ~ dstatus)
)

# Compare transversion level in WD vs other samples
with(
  patient_counts %>%
    filter(diff_status != "PD", type == "Transversion") %>%
    mutate(dstatus = ifelse(diff_status == "WD", "WD", "Other")),
  wilcox.test(prop ~ dstatus)
)

# Omnibus test for differences between WD/PD/Others
with(
  patient_counts %>%
    filter(type == "Transversion"),
  kruskal.test(prop ~ diff_status)
)

# Plot Transversion percentages per group
patient_counts %>%
  mutate(diff_status = case_when(
    diff_status == "Unknown" ~ "Other",
    TRUE ~ diff_status
  )) %>%
  mutate(diff_status = factor(diff_status, levels = c("Other", "WD", "PD"))) %>%
  filter(type == "Transversion") %>%
  ggplot(aes(diff_status, prop)) +
  geom_boxplot(aes(fill = diff_status)) +
  theme_bw() +
  labs(x = "Differentiation Level", y = "% Transversion Mutations") +
  guides(fill = FALSE)



#########################################################
## Plotting Code
#########################################################
sampleid_order <- rev(unique(patient_counts %>% filter(cohort == "South") %>% pull(Sample.ID)))

colors <- c("#3D8249", "#00A1D5")

immune_plot <- patient_counts %>%
  filter(cohort == "South") %>%
  mutate(ImmuneStatus = fct_recode(
    ImmuneStatus,
    CardiacTransplant = "CT",
    ImmunoCompetent = "IC",
    ImmunoSuppressed = "IS",
    RenalTransplant = "RT"
  )) %>%
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
  labs(x = "SCC-Normal Pairs", y = "% Mutations") +
  facet_wrap(~ImmuneStatus, scales = "free_x")
print(immune_plot)

diff_plot <- patient_counts %>%
  filter(cohort == "South") %>%
  mutate(ImmuneStatus = fct_recode(
    ImmuneStatus,
    CardiacTransplant = "CT",
    ImmunoCompetent = "IC",
    ImmunoSuppressed = "IS",
    RenalTransplant = "RT"
  )) %>%
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
  labs(x = "SCC-Normal Pairs", y = "% Mutations") +
  facet_wrap(~diff_status, scales = "free_x")
print(diff_plot)


# patient_counts %>%
#   filter(cohort == "South") %>%
#   mutate(ImmuneStatus = fct_recode(
#     ImmuneStatus,
#     CardiacTransplant = "CT",
#     ImmunoCompetent = "IC",
#     ImmunoSuppressed = "IS",
#     RenalTransplant = "RT"
#   )) %>%
#   ggplot(aes(as.factor(Sample.ID), prop, fill = type)) +
#   geom_bar(position = "fill", stat = "identity") +
#   facet_wrap(~diff_status, scales = "free_x") +
#   theme_bw() +
#   labs(x = "", y = "% Mutations") +
#   scale_fill_manual(values = colors) +
#   guides(fill = guide_legend(title = ""))


