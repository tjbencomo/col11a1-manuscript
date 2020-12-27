# File: generate_cesc_cohort.R
# Description: Clean TCGA data for CESC cohort into usable format
# for gene signature survival analysis. Use survival data
# from TCGA Curated Survival paper and phenotype information
# from Phenotypes file. All data downloaded from UCSC Xena
# December 2020.


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


cancer_type <- "CESC"
tcga_dir <- file.path("~/lee-lab-data/tcga-data", tolower(cancer_type))
col11a1_data_dir <- "data"
clinical <- read_delim(file.path(tcga_dir, 
                                 paste(cancer_type, "_clinicalMatrix", sep="")), 
                       "\t", escape_double = FALSE, trim_ws = TRUE)
survival <- read_delim(file.path(tcga_dir, paste(cancer_type, "_survival.txt.gz", sep="")), 
                       "\t", escape_double = FALSE, trim_ws = TRUE)
expr_data <- read_delim(file.path(tcga_dir, "HiSeqV2.gz"), "\t", 
                        escape_double = FALSE, trim_ws = TRUE) %>%
  column_to_rownames(var = "sample")
interested.genes <- read_excel(file.path(col11a1_data_dir, "survival_genes.xlsx"),
                               sheet = "Up") %>% pull(f)

expr_data <- reorganize_dataframe(expr_data)
symbol_conversion <- c(EOGT="C3orf64", PRRC2C="BAT2L2", CTDNEP1="DULLARD", 
                       ARHGEF28="RGNEF", DUS2="DUS2L", AKIP1="C11orf17", 
                       ATG13="KIAA0652", CTSV="CTSL2", MISP="C19orf21", 
                       SEPTIN9="SEPT9")
interested.genes[interested.genes == "43717"] <- "SEPTIN9"
interested.genes <- interested.genes[!interested.genes %in% c("RP11-231C14.4", "SMIM22")]
expr_data <- expr_data %>%
  rename(!!symbol_conversion) %>%
  select(c("Sample", all_of(interested.genes)))

patients <- expr_data %>%
  select(Sample, all_of(interested.genes)) %>%
  left_join(clinical, by=c("Sample"="sampleID")) %>%
  rename(age = "age_at_initial_pathologic_diagnosis") %>%
  filter(sample_type == "Primary Tumor") %>%
  distinct(patient_id, .keep_all = T) %>%
  select(Sample, patient_id, `_PATIENT`, age, gender, clinical_stage,
         radiation_therapy, all_of(interested.genes)) %>%
  drop_na() %>%
  mutate(clinical_stage = forcats::fct_collapse(
    clinical_stage,
    "Stage I" = c("Stage I", "Stage IA", "Stage IA1",  "Stage IA2","Stage IB", "Stage IB1", "Stage IB2"),
    "Stage II" = c("Stage II", "Stage IIA", "Stage IIA1", "Stage IIA2", "Stage IIB"),
    "Stage III" = c("Stage III", "Stage IIIA", "Stage IIIB"),
    "Stage IV" = c("Stage IVA", "Stage IVB")
  )) %>%
  inner_join(survival) %>%
  select(Sample, patient_id, `_PATIENT`, OS, OS.time, age, gender,
         clinical_stage, radiation_therapy, all_of(interested.genes))

outfp <- outfp <- file.path("data", "tcga_cesc_cohort.csv")
write_csv(patients, outfp)
