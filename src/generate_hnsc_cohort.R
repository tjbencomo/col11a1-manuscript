# File: generate_hnsc_cohort.R
# Author: Tomas Bencomo
# Email: tjbencomo@gmail.com
# Description: This file creates a finalized
# CSV file with all the clinical and transcriptomic
# info for survival analysis. The data wrangling
# code in this file is from hnsc_revised_analysis.Rmd.
# TCGA HNSC data downloaded from UCSC Xena is cleaned
# and a subset of patients selected for analysis. 
# These datasets were downloaded from the Phenotypes
# section, not the Curated Survival Data section.

library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)


reorganize_dataframe <- function(df) {
  # df should have the genes as rownames and not a column
  df_t <- data.table::transpose(df)
  colnames(df_t) <- rownames(df)
  rownames(df_t) <- colnames(df)
  df_t <- tibble::rownames_to_column(df_t, var = "Sample")
}

# Loads TCGA data from given directory and extracts relevant 
# clinical/transcriptomic data
load_data <- function(tcga_dir, local_data_dir, cancer) {
  clinical <- read_delim(file.path(tcga_dir, 
                                   paste(cancer_type, "_clinicalMatrix", sep="")), 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
  expr_data <- read_delim(file.path(tcga_dir, "HiSeqV2"), "\t", 
                           escape_double = FALSE, trim_ws = TRUE) %>%
    column_to_rownames(var = "sample")
  interested.genes <- read_excel(file.path(local_data_dir, "survival_genes.xlsx"), 
                                 sheet = "Up") %>% pull(f)
  expr_data <- reorganize_dataframe(expr_data)
  
  # Missing gene names were manually investigated to look for aliases
  # 2 genes could not be found, leaving 264 of 266 genes for the investigation
  # See hnsc_revised_analysis.Rmd for more info
  symbol_conversion <- c(EOGT="C3orf64", PRRC2C="BAT2L2", CTDNEP1="DULLARD", 
                         ARHGEF28="RGNEF", DUS2="DUS2L", AKIP1="C11orf17", 
                         ATG13="KIAA0652", CTSV="CTSL2", MISP="C19orf21", 
                         SEPTIN9="SEPT9")
  interested.genes[interested.genes == "43717"] <- "SEPTIN9"
  interested.genes <- interested.genes[!interested.genes %in% c("RP11-231C14.4", "SMIM22")]
  expr_data <- expr_data %>%
    rename(!!symbol_conversion) %>%
    select(c("Sample", interested.genes))
  
  # Filter patients based on missing data and sample type/duplicate patients
  patients <- expr_data %>%
    select(Sample, interested.genes) %>%
    left_join(clinical, by=c("Sample"="sampleID")) %>%
    mutate(age = age_at_initial_pathologic_diagnosis) %>%
    filter(sample_type == "Primary Tumor") %>%
    distinct(patient_id, .keep_all = T) %>%
    select(Sample, patient_id, `_PATIENT`, OS, OS.time, age, gender, clinical_stage,
           radiation_therapy, interested.genes) %>%
    drop_na() %>%
    mutate(clinical_stage = forcats::fct_collapse(
      clinical_stage,
      "Stage I" = c("Stage I"),
      "Stage II" = c("Stage II"),
      "Stage III" = c("Stage III"),
      "Stage IV" = c("Stage IVA", "Stage IVB", "Stage IVC")
      )
    )
}

cancer_type <- "HNSC"
tcga_dir <- file.path("~/lee-lab-data/tcga_data/", tolower(cancer_type))
data_dir <- "data"
patients <- load_data(tcga_dir, data_dir, cancer_type)

write_csv(patients, "data/tcga_hnsc_cohort.csv")
