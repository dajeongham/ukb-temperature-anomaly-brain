
###############################################################
## 02_main_analysis.R
##
## Primary association analysis: multivariable linear regression of
## longitudinal brain volume change on cumulative anomaly-day exposure.
## Sequentially adjusted models (Crude, +Demographic, +Lifestyle,
## +Scan site [primary], +Health), each additionally adjusted for the
## mean temperature anomaly. Benjamini-Hochberg FDR across the 20
## outcome-exposure combinations.
##
## Input:  data/mri_dcnt_daily.csv
## Output: output/results_model_comparison.csv
###############################################################
library(data.table)

###############################################################
## 1. Load data
###############################################################

dt <- fread("data/mri_dcnt_daily.csv")
cat("Participants:", nrow(dt), "| Variables:", ncol(dt), "\n")

###############################################################
## 2. Preprocessing
###############################################################

na_like <- c("Prefer not to answer","Do not know","None of the above","")

dt[, p21003_i2 := as.numeric(p21003_i2)]  # age
dt[, p21001_i2 := as.numeric(p21001_i2)]  # BMI
dt[, p22189    := as.numeric(p22189)]      # Townsend deprivation index
dt[, SBP       := as.numeric(p4080_i2_a0)] # systolic blood pressure
dt[, DBP       := as.numeric(p4079_i2_a0)] # diastolic blood pressure

dt[, p31 := factor(p31, levels=c("Male","Female"))]

dt[p20116_i2 %in% na_like, p20116_i2 := NA_character_]
dt[, p20116_i2 := factor(p20116_i2, levels=c("Never","Previous","Current"))]

dt[p1558_i2 %in% na_like, p1558_i2 := NA_character_]
dt[, p1558_i2 := factor(p1558_i2,
                        levels=c("Never","Special occasions only","One to three times a month",
                                 "Once or twice a week","Three or four times a week",
                                 "Daily or almost daily"))]

dt[p2443_i2 %in% na_like, p2443_i2 := NA_character_]
dt[, p2443_i2 := factor(ifelse(p2443_i2=="Yes","Yes","No"), levels=c("No","Yes"))]

# physical activity (any)
if (!"pa_any" %in% names(dt)) {
  raw_pa  <- dt$p6164_i2
  raw_pa[raw_pa %in% na_like] <- NA_character_
  pipe_pa <- ifelse(is.na(raw_pa), NA_character_, paste0("|",raw_pa,"|"))
  pa_choices <- c(
    "Walking for pleasure (not as a means of transport)",
    "Other exercises (eg: swimming, cycling, keep fit, bowling)",
    "Strenuous sports","Light DIY (eg: pruning, watering the lawn)",
    "Heavy DIY (eg: weeding, lawn mowing, carpentry, digging)")
  pa_hits <- rowSums(sapply(pa_choices, function(ch)
    !is.na(pipe_pa) & grepl(paste0("|",ch,"|"), pipe_pa, fixed=TRUE)))
  dt[, pa_any := as.integer(pa_hits > 0)]
}

dt[p21000_i0 %in% c("Prefer not to answer","Do not know"),
   p21000_i0 := NA_character_]
if (!"ethnicity" %in% names(dt)) {
  dt[, ethnicity := fcase(
    p21000_i0 %in% c("White","Any other white background","British","Irish"), "White",
    p21000_i0 %in% c("Asian or Asian British","Any other Asian background",
                     "Indian","Pakistani","Bangladeshi","Chinese"),           "Asian",
    p21000_i0 %in% c("Black or Black British","Any other Black background",
                     "African","Caribbean"),                                  "Black",
    p21000_i0 %in% c("Mixed","Any other mixed background","White and Asian",
                     "White and Black African","White and Black Caribbean"),  "Mixed",
    p21000_i0 == "Other ethnic group",                                        "Other",
    default=NA_character_)]
}
dt[, ethnicity := factor(ethnicity, levels=c("White","Asian","Black","Mixed","Other"))]

dt[p6138_i0 %in% c("Prefer not to answer","Do not know"), p6138_i0 := NA_character_]
if (!"edu_grouped" %in% names(dt)) {
  dt[, edu_grouped := fcase(
    grepl("College or University degree",     p6138_i0, fixed=TRUE), "Degree",
    grepl("A levels/AS levels or equivalent", p6138_i0, fixed=TRUE) |
      grepl("Other professional qualifications",p6138_i0, fixed=TRUE), "A-level",
    grepl("O levels/GCSEs or equivalent",    p6138_i0, fixed=TRUE) |
      grepl("CSEs or equivalent",             p6138_i0, fixed=TRUE) |
      grepl("NVQ or HND or HNC or equivalent",p6138_i0, fixed=TRUE), "O-level",
    default=NA_character_)]
}
dt[, edu_grouped := factor(edu_grouped, levels=c("O-level","A-level","Degree"))]
dt[, p54_i2 := factor(p54_i2)]
dt[, p54_i3 := factor(p54_i3)]

cat("Preprocessing done!\n")

###############################################################
## 3. Covariate definitions
###############################################################

covar_m1 <- NULL
covar_m2 <- c("p21003_i2","p31","ethnicity","p22189","edu_grouped")
covar_m3 <- c(covar_m2, "p20116_i2","p1558_i2","pa_any")
covar_m4 <- c(covar_m3, "followup_years","p54_i2","p54_i3")   # primary model
covar_m5 <- c(covar_m4, "p21001_i2","p2443_i2")               # sensitivity

###############################################################
## 4. Exposure variable definitions
###############################################################

outcomes <- c("delta_brain_norm","delta_gm_norm","delta_wm_norm",
              "delta_wmh_log","delta_hipp_bilateral")

exposures_main <- c(
  "warm_mod_days_1yr","warm_ext_days_1yr","cold_mod_days_1yr","cold_ext_days_1yr",
  "warm_mod_days_2yr","warm_ext_days_2yr","cold_mod_days_2yr","cold_ext_days_2yr",
  "warm_mod_days_3yr","warm_ext_days_3yr","cold_mod_days_3yr","cold_ext_days_3yr"
)

###############################################################
## 5. Regression function
###############################################################

run_lm <- function(outcome, exposure, dt, covs, model_name) {
  vars  <- intersect(c(outcome, exposure, covs), names(dt))
  dt_cc <- dt[complete.cases(dt[, ..vars])]
  if (nrow(dt_cc) < 100) return(NULL)

  fm <- if (is.null(covs) || length(covs)==0) {
    as.formula(paste(outcome, "~", exposure))
  } else {
    as.formula(paste(outcome, "~", exposure, "+", paste(covs, collapse=" + ")))
  }

  fit <- lm(fm, data=dt_cc)
  s   <- summary(fit)
  if (!exposure %in% rownames(s$coefficients)) return(NULL)
  cr  <- s$coefficients[exposure,]

  data.table(
    model    = model_name,
    outcome  = outcome,
    exposure = exposure,
    n        = nrow(dt_cc),
    beta     = round(cr["Estimate"],                         4),
    se       = round(cr["Std. Error"],                       4),
    p_val    = round(cr["Pr(>|t|)"],                         4),
    ci_lo    = round(cr["Estimate"] - 1.96*cr["Std. Error"], 4),
    ci_hi    = round(cr["Estimate"] + 1.96*cr["Std. Error"], 4)
  )
}

run_all <- function(exposures, covs, model_name) {
  cat("\n[", model_name, "] running...\n")
  rbindlist(lapply(outcomes, function(out)
    rbindlist(lapply(exposures, function(exp)
      tryCatch(run_lm(out, exp, dt, covs, model_name),
               error=function(e) NULL)
    ), fill=TRUE)
  ), fill=TRUE)
}

# Anomaly-adjusted version (adds the mean_dcnt covariate)
run_all_dcnt_adj <- function(exposures, covs, model_name) {
  cat("\n[", model_name, "] running...\n")
  rbindlist(lapply(outcomes, function(out)
    rbindlist(lapply(exposures, function(exp) {
      win      <- regmatches(exp, regexpr("\\dyr$", exp))
      covs_adj <- c(covs, paste0("mean_dcnt_", win))
      tryCatch(run_lm(out, exp, dt, covs_adj, model_name),
               error=function(e) NULL)
    }), fill=TRUE)
  ), fill=TRUE)
}

###############################################################
## 6. Phase 1: run M1-M5
###############################################################

cat("\n===== Phase 1: M1-M5 model comparison =====\n")

results_model <- rbindlist(list(
  run_all(exposures_main, covar_m1, "M1_crude"),
  run_all(exposures_main, covar_m2, "M2_demographic"),
  run_all(exposures_main, covar_m3, "M3_lifestyle"),
  run_all(exposures_main, covar_m4, "M4_scan"),
  run_all(exposures_main, covar_m5, "M5_full")
), fill=TRUE)

results_model[, p_fdr := p.adjust(p_val, method="BH"), by=.(outcome, model)]
results_model[, window := sub(".*_days_", "", exposure)]
results_model[, p_fdr_20 := p.adjust(p_val, method = "BH"),
                 by = .(model, window)]

###############################################################
## 7. Phase 2: anomaly-adjusted version (M1b-M5b)
###############################################################

cat("\n===== Phase 2: anomaly-adjusted version =====\n")

results_dcnt_adj <- rbindlist(list(
  run_all_dcnt_adj(exposures_main, covar_m1, "M1b_crude_dcntadj"),
  run_all_dcnt_adj(exposures_main, covar_m2, "M2b_demo_dcntadj"),
  run_all_dcnt_adj(exposures_main, covar_m3, "M3b_life_dcntadj"),
  run_all_dcnt_adj(exposures_main, covar_m4, "M4b_scan_dcntadj"),   # primary
  run_all_dcnt_adj(exposures_main, covar_m5, "M5b_full_dcntadj")    # sensitivity
), fill=TRUE)

results_dcnt_adj[, p_fdr := p.adjust(p_val, method="BH"), by=.(outcome, model)]
results_dcnt_adj[, window := sub(".*_days_", "", exposure)]
results_dcnt_adj[, p_fdr_20 := p.adjust(p_val, method = "BH"),
            by = .(model, window)]

###############################################################
## 8. Combine & save
###############################################################

results_all <- rbindlist(list(results_model, results_dcnt_adj), fill=TRUE)

cat("\n===== Results summary =====\n")
cat("\n--- M4 (primary) FDR q<0.05 ---\n")
print(results_all[model=="M4_scan" & p_fdr<0.05,
                  .(outcome, exposure, n, beta, p_val, p_fdr)][order(p_fdr)])
cat("\n--- M4b (anomaly-adjusted, primary) FDR q<0.05 ---\n")
print(results_all[model=="M4b_scan_dcntadj" & p_fdr<0.05,
                  .(outcome, exposure, n, beta, p_val, p_fdr)][order(p_fdr)])

out_dir <- "output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
fwrite(results_all, file.path(out_dir, "results_model_comparison.csv"))
cat("\nSaved: results_model_comparison.csv\n")
