###############################################################
## 02_main_analysis.R
##
## Primary association analysis: multivariable linear regression of
## longitudinal brain volume change on cumulative anomaly-day exposure.
##
## Five sequentially adjusted models:
##   Crude        - exposure only
##   +Demographic - age, sex, ethnicity, Townsend deprivation
##   +Lifestyle   - + education, smoking, alcohol, physical activity
##   +Scan site   - + imaging centre, follow-up duration   [PRIMARY]
##   +Health      - + BMI, diabetes                         [sensitivity]
## All models except Crude additionally adjust for the corresponding
## window's mean temperature anomaly (continuous), to separate warm/cold
## anomalies from broader warming/cooling trends.
##
## Each exposure category is fitted in a separate model (the four nested
## categories are never entered jointly). Benjamini-Hochberg FDR is applied
## across the 20 outcome-exposure combinations (5 outcomes x 4 categories)
## within each exposure window.
##
## Input:  data/mri_dcnt_daily.csv
## Output: output/results_model_comparison.csv
###############################################################
library(data.table)

dt <- fread("data/mri_dcnt_daily.csv")

###############################################################
## 1. Preprocessing
###############################################################

na_like <- c("Prefer not to answer","Do not know","None of the above","")

dt[, p21003_i2 := as.numeric(p21003_i2)]   # age
dt[, p21001_i2 := as.numeric(p21001_i2)]   # BMI
dt[, p22189    := as.numeric(p22189)]      # Townsend deprivation index
dt[, p31       := factor(p31, levels=c("Male","Female"))]

dt[p20116_i2 %in% na_like, p20116_i2 := NA_character_]
dt[, p20116_i2 := factor(p20116_i2, levels=c("Never","Previous","Current"))]

dt[p1558_i2 %in% na_like, p1558_i2 := NA_character_]
dt[, p1558_i2 := factor(p1558_i2,
     levels=c("Never","Special occasions only","One to three times a month",
              "Once or twice a week","Three or four times a week",
              "Daily or almost daily"))]

dt[p2443_i2 %in% na_like, p2443_i2 := NA_character_]
dt[, p2443_i2 := factor(ifelse(p2443_i2=="Yes","Yes","No"), levels=c("No","Yes"))]

# physical activity (any of the listed activity types)
raw_pa  <- dt$p6164_i2; raw_pa[raw_pa %in% na_like] <- NA_character_
pa_choices <- c(
  "Walking for pleasure (not as a means of transport)",
  "Other exercises (eg: swimming, cycling, keep fit, bowling)",
  "Strenuous sports","Light DIY (eg: pruning, watering the lawn)",
  "Heavy DIY (eg: weeding, lawn mowing, carpentry, digging)")
dt[, pa_any := as.integer(rowSums(sapply(pa_choices, function(ch)
     !is.na(raw_pa) & grepl(ch, raw_pa, fixed=TRUE))) > 0)]

dt[p21000_i0 %in% c("Prefer not to answer","Do not know"), p21000_i0 := NA_character_]
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
dt[, ethnicity := factor(ethnicity, levels=c("White","Asian","Black","Mixed","Other"))]

dt[p6138_i0 %in% c("Prefer not to answer","Do not know"), p6138_i0 := NA_character_]
dt[, edu_grouped := fcase(
  grepl("College or University degree", p6138_i0, fixed=TRUE), "Degree",
  grepl("A levels/AS levels or equivalent", p6138_i0, fixed=TRUE) |
    grepl("Other professional qualifications", p6138_i0, fixed=TRUE), "A-level",
  grepl("O levels/GCSEs or equivalent", p6138_i0, fixed=TRUE) |
    grepl("CSEs or equivalent", p6138_i0, fixed=TRUE) |
    grepl("NVQ or HND or HNC or equivalent", p6138_i0, fixed=TRUE), "O-level",
  default=NA_character_)]
dt[, edu_grouped := factor(edu_grouped, levels=c("O-level","A-level","Degree"))]

dt[, p54_i2 := factor(p54_i2)]   # imaging centre (t0)
dt[, p54_i3 := factor(p54_i3)]   # imaging centre (t1)

###############################################################
## 2. Outcomes, exposures, and sequential covariate sets
###############################################################

outcomes <- c("delta_brain_norm","delta_gm_norm","delta_wm_norm",
              "delta_wmh_log","delta_hipp_bilateral")

# 4 nested categories x 3 windows (1-, 2-, 3-year)
exposures <- as.vector(outer(
  c("warm_mod_days","warm_ext_days","cold_mod_days","cold_ext_days"),
  c("1yr","2yr","3yr"), paste, sep="_"))

demo <- c("p21003_i2","p31","ethnicity","p22189")
life <- c(demo, "edu_grouped","p20116_i2","p1558_i2","pa_any")
scan <- c(life, "followup_years","p54_i2","p54_i3")
hlth <- c(scan, "p21001_i2","p2443_i2")

# anom = TRUE adds the window-specific mean temperature anomaly
models <- list(
  Crude        = list(covs = NULL, anom = FALSE),
  Demographic  = list(covs = demo, anom = TRUE),
  Lifestyle    = list(covs = life, anom = TRUE),
  ScanSite     = list(covs = scan, anom = TRUE),   # primary
  Health       = list(covs = hlth, anom = TRUE)    # sensitivity
)

###############################################################
## 3. Fit one outcome x exposure x model
###############################################################

fit_one <- function(outcome, exposure, model_name) {
  m    <- models[[model_name]]
  win  <- sub(".*_days_", "", exposure)                 # "1yr"/"2yr"/"3yr"
  covs <- m$covs
  if (m$anom) covs <- c(covs, paste0("mean_dcnt_", win))

  vars <- intersect(c(outcome, exposure, covs), names(dt))
  d    <- dt[complete.cases(dt[, ..vars])]
  if (nrow(d) < 100) return(NULL)

  cf <- summary(lm(reformulate(c(exposure, covs), outcome), data=d))$coefficients
  if (!exposure %in% rownames(cf)) return(NULL)
  b <- cf[exposure, ]

  data.table(model=model_name, window=win, outcome=outcome, exposure=exposure,
             n=nrow(d), beta=round(b[1],4), se=round(b[2],4),
             ci_lo=round(b[1]-1.96*b[2],4), ci_hi=round(b[1]+1.96*b[2],4),
             p=round(b[4],4))
}

###############################################################
## 4. Run all combinations and apply FDR within each window
###############################################################

grid <- CJ(outcome=outcomes, exposure=exposures, model=names(models), sorted=FALSE)
results <- rbindlist(Map(fit_one, grid$outcome, grid$exposure, grid$model))

# BH-FDR across the 20 outcome-exposure combinations within each (model, window)
results[, q := p.adjust(p, method="BH"), by=.(model, window)]

###############################################################
## 5. Save and report the primary model
###############################################################

dir.create("output", showWarnings=FALSE, recursive=TRUE)
fwrite(results, "output/results_model_comparison.csv")

cat("Primary model (+Scan site), 1-year window, q < 0.05:\n")
print(results[model=="ScanSite" & window=="1yr" & q < 0.05,
              .(outcome, exposure, n, beta, p, q)][order(q)])
cat("\nSaved: output/results_model_comparison.csv\n")
