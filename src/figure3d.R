# Description: Plot CESC Kaplan Meier plot with Mut COL11A1 Signature


library(readr)
library(dplyr)
library(ggplot2)
library(rms)
library(survminer)

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

ddist <- datadist(df)
options(datadist = "ddist")

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


### Plotting Code for Manuscript
km_df <- df %>%
  mutate(cut_sig = ifelse(signature > 0, "Upregulated", "Downregulated")) %>%
  mutate(cut_sig = factor(cut_sig, levels = c("Upregulated", "Downregulated")))
km.fit <- survfit(Surv(OS.time, OS) ~ cut_sig, data=km_df)
km_plot <- ggsurvplot(km.fit, data=km_df, legend.lab = c("High Expression", "Low Expression"))$plot + 
  labs(x = "Time (Days)", y = "CESC (TCGA) Survival Probability (%)")
print(km_plot)
