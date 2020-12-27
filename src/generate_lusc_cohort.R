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

load_data <- function(tcga_dir, local_data_dir, cancer) {
  clinical <- read_delim(file.path(tcga_dir, 
                                   paste(cancer_type, "_clinicalMatrix", sep="")), 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
  expr_data <- read_delim(file.path(tcga_dir, "HiSeqV2.gz"), "\t", 
                          escape_double = FALSE, trim_ws = TRUE) %>%
    column_to_rownames(var = "sample")
  # interested.genes <- read_excel(file.path(local_data_dir, "survival_genes.xlsx"), 
  #                                sheet = "Up") %>% pull(f)
  interested.genes <- scan(file.path(local_data_dir, "264_gene_signature.txt"),
                           what = character())
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
    select(Sample, all_of(interested.genes)) %>%
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

# cancer_type <- "LUSC"
# tcga_dir <- file.path("~/lee-lab-data/tcga-data", tolower(cancer_type))
# data_dir <- "data"
# patients <- load_data(tcga_dir, data_dir, cancer_type)

cancer_type <- "LUSC"
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
  select(Sample, patient_id, `_PATIENT`, age, gender, pathologic_stage,
         radiation_therapy, all_of(interested.genes)) %>%
  drop_na() %>%
  mutate(pathologic_stage = forcats::fct_collapse(
    pathologic_stage,
    "Stage I" = c("Stage I", "Stage IA", "Stage IB"),
    "Stage II" = c("Stage II", "Stage IIA", "Stage IIB"),
    "Stage III" = c("Stage III", "Stage IIIA", "Stage IIIB"),
    "Stage IV" = c("Stage IV"),
    "discrepancy" = c("[Discrepancy]")
  )) %>%
  inner_join(survival) %>%
  select(Sample, patient_id, `_PATIENT`, OS, OS.time, age, gender,
         pathologic_stage, radiation_therapy, all_of(interested.genes))

outfp <- file.path("data", "tcga_lusc_cohort.csv")
write_csv(patients, outfp)


