# File: cesc_survival_analysis.R
# Description: Investigate correlation between
# 264 gene signature from COL11A1 RNA-Seq experiments
# and overall survival in CESC. 


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
ggsurvplot(km.fit, data=km_df, legend.lab = c("High Expression", "Low Expression"))$plot + 
  labs(x = "Time (Days)", y = "CESC (TCGA) Survival Probability (%)")


## Marginal Effect Plot
full.rms.fit <- cph(Surv(OS.time, OS) ~ age + clinical_stage + 
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
table_fp <- file.path(data_dir, "cesc_cox_table.csv")
write_csv(fit.table, table_fp)



