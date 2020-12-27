# File: cesc_survival_analysis.R
# Description: Investigate correlation between
# 264 gene signature from COL11A1 RNA-Seq experiments
# and overall survival in CESC. 


library(readr)
library(dplyr)
library(ggplot2)
library(rms)

data_dir <- "data"
patient_fp <- file.path(data_dir, "tcga_cesc_cohort.csv")
genesig_fp <- file.path(data_dir, "264_gene_signature.txt")
patients <- read_csv(patient_fp)
genes <- scan(genesig_fp, what = character())

scaled_expr <- scale(patients[, genes])

df <- cbind(
  select(patients, OS.time, OS, age, gender, clinical_stage,
         radiation_therapy),
  scaled_expr
) %>%
  mutate(
    signature = scale(rowSums(.[genes]))[, 1]
  )

sig.fit <- coxph(Surv(OS.time, OS) ~ signature, data=df)
print("264 Gene signature assocation with overall survival:")
print(summary(sig.fit))
print("Proportional hazards check")
print(cox.zph(sig.fit))

full.fit <- coxph(Surv(OS.time, OS) ~ age + clinical_stage + 
                    strata(radiation_therapy) + signature, data=df)
print("Association between gene signature and OS after adjustments")
print(summary(full.fit))
print("Proportional hazards check")
print(cox.zph(full.fit))
