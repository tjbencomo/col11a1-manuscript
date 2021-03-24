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
library(survminer)

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
  ) %>%
  filter(pathologic_stage != "discrepancy")

ddist <- datadist(df)
options(datadist = "ddist")

sig.fit <- coxph(Surv(OS.time, OS) ~ signature, data=df)
print("264 Gene signature assocation with overall survival:")
print(summary(sig.fit))
print("Proportional hazards check")
print(cox.zph(sig.fit))

sig.rcs.fit <- cph(Surv(OS.time, OS) ~ rcs(signature, 3), data=df,
                   x=T, y=T)
print(anova(sig.rcs.fit))
print(cox.zph(sig.rcs.fit))
ggplot(Predict(sig.rcs.fit)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "264 Gene Signature (Z-Score)",
       y = "Log Hazard Ratio",
       caption = "")

full.fit <- coxph(Surv(OS.time, OS) ~ age + gender + pathologic_stage + 
                    strata(radiation_therapy) + signature, data=df,
                  x=T, y=T)
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

## Likelihood Ratio Test and Added Info
subset.fit <- cph(Surv(OS.time, OS) ~ age + gender + pathologic_stage
                  + strat(radiation_therapy), data=df,
                  x=T, y=T)
full.rms.fit <- cph(Surv(OS.time, OS) ~ age + gender + pathologic_stage + 
                        strat(radiation_therapy) + signature, data=df,
                      x=T, y=T)
print(lrtest(subset.fit, full.rms.fit))
print(paste("New information added: ", 
            1 - (subset.fit$stats['Model L.R.'] / full.rms.fit$stats['Model L.R.']),
            "%", sep=""))

## Plotting Code for Manuscript

## KM Plot
km_df <- df %>%
  mutate(cut_sig = ifelse(signature > 0, "Upregulated", "Downregulated")) %>%
  mutate(cut_sig = factor(cut_sig, levels = c("Upregulated", "Downregulated")))
km.fit <- survfit(Surv(OS.time, OS) ~ cut_sig, data=km_df)
ggsurvplot(km.fit, data=km_df, legend.lab = c("High Expression", "Low Expression"))$plot +
  labs(x = "Time (Days)", y = "LUSC (TCGA) Survival Probability (%)")


## Marginal Effect Plot
full.rms.fit <- cph(Surv(OS.time, OS) ~ age + gender + pathologic_stage + 
                    strat(radiation_therapy) + signature, data=df,
                  x=T, y=T)
ggplot(Predict(full.rms.fit, signature)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "264 Gene Signature (Z-Score)",
       y = "Log Hazard Ratio",
       caption = "")

## Cox Regression Table
fit.table <- broom::tidy(full.fit) %>%
  mutate(HR = exp(estimate),
         low.ci = exp(conf.low),
         high.ci = exp(conf.high),
         CI = paste(round(low.ci, 2), "-", round(high.ci, 2), sep = "")) %>%
  select(term, HR, CI, p.value) %>%
  rename(`Prognostic Factor` = "term", `P-Value` = p.value)
fit.table
table_fp <- file.path(data_dir, "lusc_cox_table.csv")
write_csv(fit.table, table_fp)
