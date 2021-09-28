library(readr)
library(dplyr)
library(rms)
library(survminer)
library(broom)
library(ggplot2)

hnsc_analysis <- function() {
  print("Analyzing HNSC Cohort")
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
  print(ddist)
  
  full_sig.fit <- cph(Surv(OS.time, OS) ~ age + gender + clinical_stage + 
                        radiation_therapy + full_signature, data=patients, x=T, y=T)
  print(full_sig.fit)
  
  full_sig.phtest <- cox.zph(full_sig.fit)
  print(full_sig.phtest)
  
  full_sig.fit.strat <- cph(Surv(OS.time, OS) ~ age + gender + clinical_stage + 
                              strat(radiation_therapy) + full_signature, data=patients, x=T, y=T)
  print(full_sig.fit.strat)
  
  print(anova(full_sig.fit.strat))
  
  full_sig.phtest.strat <- cox.zph(full_sig.fit.strat)
  print(full_sig.phtest.strat)
  
  subset.fit <- cph(Surv(OS.time, OS) ~ age + gender + clinical_stage + 
                      strat(radiation_therapy), data=patients, x=T, y=T)
  print(lrtest(subset.fit, full_sig.fit.strat))
  print(paste("New information added: ", 
              (1 - (subset.fit$stats['Model L.R.'] / full_sig.fit.strat$stats['Model L.R.'])) * 100,
              "%", sep=""))
  
  effect_plot <- ggplot(Predict(full_sig.fit.strat, 
                 full_signature)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ggtitle("Association Between 264 Gene Signature And Overall Survival") +
    labs(x = "Pathway Expression (Z-Score)",
         y = "Log Hazard Ratio",
         caption = "") +
    annotate("text", x = -1.5, y = .75, label = "p = 5.56e-4")
  
  full_sig.coxph.fit <- coxph(Surv(OS.time, OS) ~ age + gender + clinical_stage + 
                                strata(radiation_therapy) + full_signature, data=patients)
  full_sig.fit.table <- broom::tidy(full_sig.coxph.fit) %>%
    mutate(HR = exp(estimate),
           low.ci = exp(estimate - (1.96 *std.error)),
           high.ci = exp(estimate + (1.96 *std.error)),
           CI = paste(round(low.ci, 2), "-", round(high.ci, 2), sep = "")) %>%
    select(term, HR, CI, p.value) %>%
    rename(`Prognostic Factor` = "term", `P-Value` = p.value)
  return(list("plot" = effect_plot, "table" = full_sig.fit.table))
}

cesc_analysis <- function() {
  print("Analyzing CESC Cohort")
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
  options("datadist" = ddist)
  
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
  
  ## Likelihood Ratio Test and Added Info
  subset.fit <- cph(Surv(OS.time, OS) ~ age + clinical_stage + 
                      strat(radiation_therapy), data=df, x=T, y=T)
  full.rms.fit <- cph(Surv(OS.time, OS) ~ age + clinical_stage + 
                        strat(radiation_therapy) + signature, data=df, x=T, y=T)
  print(lrtest(subset.fit, full.rms.fit))
  print(paste("New information added: ", 
              (1 - (subset.fit$stats['Model L.R.'] / full.rms.fit$stats['Model L.R.'])) * 100,
              "%", sep=""))
  
  ## Marginal Effect Plot
  full.rms.fit <- cph(Surv(OS.time, OS) ~ age + clinical_stage + 
                        strat(radiation_therapy) + signature, data=df,
                      x=T, y=T)
  effect_plot <- ggplot(Predict(full.rms.fit, signature)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(x = "264 Gene Signature (Z-Score)",
         y = "Log Hazard Ratio",
         caption = "") 
  
  ## Cox Regression Table
  fit.table <- broom::tidy(full.fit) %>%
    mutate(HR = exp(estimate),
           low.ci = exp(estimate - (1.96 *std.error)),
           high.ci = exp(estimate + (1.96 *std.error)),
           CI = paste(round(low.ci, 2), "-", round(high.ci, 2), sep = "")) %>%
    select(term, HR, CI, p.value) %>%
    rename(`Prognostic Factor` = "term", `P-Value` = p.value)
  return(list("plot" = effect_plot, "table" = fit.table))
}

lusc_analysis <- function() {
  # Investigate correlation between
  # 264 gene signature identified from COL11A1 RNA-Seq
  # experiments and overall survival in LUSC. Unlike
  # HNSC and CESC, the gene signature  shows weak
  # evidence that suggests the PH assumption is not
  # satisfied. Restricted cubic spline modeling of the
  # gene signature corrects the PH assumption and reaches
  # the same overall conclusion, so we keep the gene signature
  # as a linear relationship to fit with the HNSC and CESC analyses.

  print("Analyzing LUSC Cohort")  
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
  options("datadist" = ddist)
  
  sig.fit <- coxph(Surv(OS.time, OS) ~ signature, data=df)
  print("264 Gene signature assocation with overall survival:")
  print(summary(sig.fit))
  print("Proportional hazards check")
  print(cox.zph(sig.fit))
  
  sig.rcs.fit <- cph(Surv(OS.time, OS) ~ rcs(signature, 3), data=df,
                     x=T, y=T)
  print(anova(sig.rcs.fit))
  print(cox.zph(sig.rcs.fit))
  # ggplot(Predict(sig.rcs.fit)) +
  #   geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  #   labs(x = "264 Gene Signature (Z-Score)",
  #        y = "Log Hazard Ratio",
  #        caption = "")
  
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
  # ggplot(Predict(full.rcs.fit, signature)) +
  #   geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  #   labs(x = "264 Gene Signature (Z-Score)",
  #        y = "Log Hazard Ratio",
  #        caption = "") +
  #   ggtitle("Gene Signature vs OS after adjustments")
  
  ## Likelihood Ratio Test and Added Info
  subset.fit <- cph(Surv(OS.time, OS) ~ age + gender + pathologic_stage
                    + strat(radiation_therapy), data=df,
                    x=T, y=T)
  full.rms.fit <- cph(Surv(OS.time, OS) ~ age + gender + pathologic_stage + 
                        strat(radiation_therapy) + signature, data=df,
                      x=T, y=T)
  print(lrtest(subset.fit, full.rms.fit))
  print(paste("New information added: ", 
              (1 - (subset.fit$stats['Model L.R.'] / full.rms.fit$stats['Model L.R.'])) * 100,
              "%", sep=""))
  
  ## Marginal Effect Plot
  full.rms.fit <- cph(Surv(OS.time, OS) ~ age + gender + pathologic_stage + 
                        strat(radiation_therapy) + signature, data=df,
                      x=T, y=T)
  effect_plot <- ggplot(Predict(full.rms.fit, signature)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(x = "264 Gene Signature (Z-Score)",
         y = "Log Hazard Ratio",
         caption = "")
  
  ## Cox Regression Table
  fit.table <- broom::tidy(full.fit) %>%
    mutate(HR = exp(estimate),
           low.ci = exp(estimate - (1.96 *std.error)),
           high.ci = exp(estimate + (1.96 *std.error)),
           CI = paste(round(low.ci, 2), "-", round(high.ci, 2), sep = "")) %>%
    select(term, HR, CI, p.value) %>%
    rename(`Prognostic Factor` = "term", `P-Value` = p.value)
  return(list("plot" = effect_plot, "table" = fit.table))
}


hnsc_results <- hnsc_analysis()
print(hnsc_results$plot)
print(hnsc_results$table)

cesc_results <- cesc_analysis()
print(cesc_results$plot)
print(cesc_results$table)

lusc_results <- lusc_analysis()
print(lusc_results$plot)
print(lusc_results$table)