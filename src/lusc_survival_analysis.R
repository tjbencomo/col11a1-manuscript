# File : lusc_survival_analysis.R
# Description: Investigate correlation between
# 264 gene signature identified from COL11A1 RNA-Seq
# experiments and overall survival in LUSC. Unlike
# HNSC and CESC, the gene signature  shows weak
# evidence that suggests the PH assumption is not
# satisfied. Restricted cubic spline modeling of the
# gene signature corrects the PH assumption and reaches
# the same overall conclusion, so we keep the gene signature
# as a linear relationship to fit with the HNSC and CESC analyses.

library(readr)
library(dplyr)
library(ggplot2)
library(rms)

data_dir <- "data"
patient_fp <- file.path(data_dir, "tcga_lusc_cohort.csv")
genesig_fp <- file.path(data_dir, "264_gene_signature.txt")
patients <- read_csv(patient_fp)
genes <- scan(genesig_fp, what = character())

scaled_expr <- scale(patients[, genes])

df <- cbind(
  select(patients, OS.time, OS, age, gender, pathologic_stage,
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

sig.rcs.fit <- cph(Surv(OS.time, OS) ~ rcs(signature, 3), data=df,
                   x=T, y=T)
anova(sig.rcs.fit)
cox.zph(sig.rcs.fit)
ggplot(Predict(sig.rcs.fit)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "264 Gene Signature (Z-Score)",
       y = "Log Hazard Ratio",
       caption = "")

full.fit <- coxph(Surv(OS.time, OS) ~ age + gender + pathologic_stage + 
                    strata(radiation_therapy) + signature, data=df)
print("Association between gene signature and OS after adjustments")
print(summary(full.fit))
print("Proportional hazards check")
print(cox.zph(full.fit))

full.rcs.fit <- cph(Surv(OS.time, OS) ~ age + gender + strat(pathologic_stage)
                    + strat(radiation_therapy) + rcs(signature, 3), data=df,
                    x=T, y=T)
anova(full.rcs.fit)
cox.zph(full.rcs.fit)
ggplot(Predict(full.rcs.fit, signature)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "264 Gene Signature (Z-Score)",
       y = "Log Hazard Ratio",
       caption = "") +
  ggtitle("Gene Signature vs OS after adjustments")
