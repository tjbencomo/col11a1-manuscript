# Description: Plot HNSC Kaplan Meier plot with Mut COL11A1 Signature

library(readr)
library(dplyr)
library(rms)
library(survminer)
library(broom)

data_dir <- "data"
patients <- read_csv(file.path(data_dir, "tcga_hnsc_cohort.csv"))

clinical_cols <- c("Sample", "patient_id", "_PATIENT", "OS", "OS.time", 
                   "age", "gender", "clinical_stage", "radiation_therapy")
full_signature <- colnames(patients)[!colnames(patients) %in% clinical_cols]

scaled_expr <- scale(patients[, full_signature])
patients <- cbind(
  select(patients, OS.time, OS, age, gender, clinical_stage,
         radiation_therapy),
  scaled_expr
) %>%
  mutate(full_signature = scale(rowSums(.[full_signature]))[, 1])

ddist <- datadist(patients)
options("datadist" = ddist)

full_sig.fit <- cph(Surv(OS.time, OS) ~ age + gender + clinical_stage + 
                      radiation_therapy + full_signature, data=patients, x=T, y=T)
full_sig.fit

full_sig.phtest <- cox.zph(full_sig.fit)
full_sig.phtest

full_sig.fit.strat <- cph(Surv(OS.time, OS) ~ age + gender + clinical_stage + 
                            strat(radiation_therapy) + full_signature, data=patients, x=T, y=T)
full_sig.fit.strat

anova(full_sig.fit.strat)

full_sig.phtest.strat <- cox.zph(full_sig.fit.strat)
full_sig.phtest.strat

patients.full_sig.km <- patients %>%
  mutate(full_sig_group = case_when(
    full_signature > 0 ~ "Upregulated",
    full_signature <= 0 ~ "Downregulated"
  )) %>%
  mutate(full_sig_group = factor(full_sig_group, levels = c("Upregulated", "Downregulated")))
full.sig.km.fit <- survfit(Surv(OS.time, OS) ~ full_sig_group, data=patients.full_sig.km)
km_plot <- ggsurvplot(full.sig.km.fit, data=patients.full_sig.km, legend.lab = c("High Expression", "Low Expression"))$plot + labs(x = "Time (Days)", y = "HNSC (TCGA) Survival Probability (%)")
print(km_plot)