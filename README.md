# ukb-temperature-anomaly-brain

Analysis code for the study:

> **Cumulative Exposure to Warm and Cold Temperature Anomalies and Longitudinal Changes in Brain Structure: A Prospective Cohort Study**

This repository contains the core R code to (1) construct the exposure —
cumulative warm and cold temperature **anomaly days**, defined as daily
departures from each location's 1961–1990 day-of-year climatological norm
exceeding ±1 and ±2 SD — and (2) estimate its primary associations with
longitudinal change in five brain structural measures (total brain, grey
matter, white matter, white matter hyperintensities [WMH], and hippocampal
volume) in the UK Biobank MRI sub-cohort.

The code is intentionally limited to the **exposure construction, data
preparation, and primary (main) regression analysis**. Secondary analyses
(subgroup, sensitivity, dose–response/GAM), tables, and figures are not
included in this public release.

## Repository structure

```
ukb-temperature-anomaly-brain/
├── R/
│   ├── 01_data_preparation.R   # exposure construction + outcomes + covariates
│   └── 02_main_analysis.R      # primary regression models + FDR
├── data/      # input data (not included; see Data availability)
├── output/    # results written here
└── README.md
```

## Data availability

Individual-level data cannot be shared. UK Biobank data are available to
approved researchers through the UK Biobank Access Management System
(https://www.ukbiobank.ac.uk). Daily ERA5-Land reanalysis temperature data are
publicly available from the Copernicus Climate Data Store. Scripts expect inputs
under `data/`; these are not distributed with this repository.

## Requirements

- R (≥ 4.3)
- R packages: `data.table`, `terra`, `ncdf4`, `lubridate`

```r
install.packages(c("data.table","terra","ncdf4","lubridate"))
```

## Pipeline

Run from the repository root, in order. Inputs are read from `data/`; results
are written to `output/`.

| Script | Purpose | Main output |
|---|---|---|
| `R/01_data_preparation.R` | Build the analysis dataset: ERA5-Land day-of-year baseline (μ_id) and cell-specific anomaly SD (σ_i); daily anomaly and classification into warm/cold (>±1 SD) and extreme (>±2 SD) anomaly days; annual and cumulative (1–3 yr) anomaly-day counts; brain volume change (t1 − t0) and covariates. | `data/mri_dcnt_daily.csv` |
| `R/02_main_analysis.R` | Primary multivariable linear regression of brain volume change on cumulative anomaly-day exposure; sequentially adjusted models (Crude, +Demographic, +Lifestyle, +Scan site [primary], +Health), each adjusted for the mean temperature anomaly; Benjamini–Hochberg FDR across the 20 outcome–exposure combinations. | `output/results_model_comparison.csv` |

## Exposure definition (Supplementary Methods S1)

For each ERA5-Land grid cell, a day-of-year baseline mean temperature (μ_id,
1961–1990) is computed, together with a cell-specific standard deviation (σ_i)
of daily anomalies over the reference period. The daily anomaly is the observed
daily mean temperature minus μ_id, and each day is classified by the magnitude
of the anomaly relative to σ_i:

- Extreme warm anomaly day: anomaly > +2 σ_i
- Warm anomaly day: anomaly > +1 σ_i
- Cold anomaly day: anomaly < −1 σ_i
- Extreme cold anomaly day: anomaly < −2 σ_i

The exposure is the cumulative number of such days in the year preceding the
baseline imaging visit (primary window); 2- and 3-year windows are also derived.

## Notes

- Some in-line comments are in Korean.
- Multiple-comparison correction uses the Benjamini–Hochberg false discovery
  rate (FDR) across the 20 outcome–exposure combinations.

## License

Code is released under the MIT License (see `LICENSE`). This applies to the
code only and grants no rights to UK Biobank or ERA5-Land data.

## Citation

If you use this code, please cite the associated publication (to be added upon
publication).
