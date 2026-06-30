
###############################################################
## 01_data_preparation.R
##
## Build the participant-level analysis dataset.
##  - ERA5-Land 1961-1990 day-of-year baseline mean and cell-specific
##    SD of daily anomalies
##  - daily anomaly = observed daily mean - baseline; classified into
##    warm/cold (>짹1 SD) and extreme warm/cold (>짹2 SD) anomaly days
##  - annual and cumulative (1-3 yr) anomaly-day counts per participant
##  - longitudinal brain volume change (t1 - t0)
##
## Output: data/mri_dcnt_daily.csv
###############################################################
library(terra)
library(ncdf4)
library(data.table)
library(lubridate)

## Paths (relative to the repository root)
mri_path        <- "data/data_brain2.csv"
era5_sample     <- "data/ERA5/t2m_2014_01.nc"          # any monthly file (for the grid)
nc_dir_baseline <- "data/ERA5/baseline_1961_1990"      # reference-period NetCDF
nc_dir_obs      <- "data/ERA5/observation"             # observation-period NetCDF
out_path        <- "data/mri_dcnt_daily.csv"

md_levels <- format(seq(as.Date("2001-01-01"), as.Date("2001-12-31"), by="day"), "%m-%d")

###############################################################
## STEP 1. Brain volume change (t1 - t0)
###############################################################

mri <- fread(mri_path)
mri[, t0_date := as.Date(p53_i2)]
mri[, t1_date := as.Date(p53_i3)]
mri[, tmp_id  := .I]

mri <- mri[!is.na(t0_date) & !is.na(t1_date) &
           !is.na(p25009_i2) & !is.na(p25009_i3) &
           !is.na(p32223_a0) & !is.na(p32224_a0)]

mri[, t0_year        := year(t0_date)]
mri[, followup_years := as.numeric(t1_date - t0_date) / 365.25]

# Five outcomes (WMH on the log scale; hippocampus = left + right)
mri[, delta_brain_norm := as.numeric(p25009_i3) - as.numeric(p25009_i2)]
mri[, delta_gm_norm    := as.numeric(p25005_i3) - as.numeric(p25005_i2)]
mri[, delta_wm_norm    := as.numeric(p25007_i3) - as.numeric(p25007_i2)]
mri[, delta_wmh_log    := log(as.numeric(p25781_i3) + 1) - log(as.numeric(p25781_i2) + 1)]
mri[, delta_hipp_bilateral :=
      (as.numeric(p25019_i3) + as.numeric(p25020_i3)) -
      (as.numeric(p25019_i2) + as.numeric(p25020_i2))]

###############################################################
## STEP 2. Residential coordinates at t0 (accounting for relocation)
###############################################################

mri[, best_array := 0L]
for (k in 1:4)
  mri[!is.na(get(paste0("p32220_a",k))) &
      as.Date(get(paste0("p32220_a",k))) <= t0_date, best_array := k]

mri[, easting := as.integer(mapply(function(a, ...) c(...)[[a+1]],
      best_array, p32223_a0, p32223_a1, p32223_a2, p32223_a3, p32223_a4))]
mri[, northing := as.integer(mapply(function(a, ...) c(...)[[a+1]],
      best_array, p32224_a0, p32224_a1, p32224_a2, p32224_a3, p32224_a4))]
mri <- mri[!is.na(easting) & !is.na(northing)]

###############################################################
## STEP 3. Match each participant to an ERA5-Land grid cell
###############################################################

nc  <- nc_open(era5_sample)
lon <- ncvar_get(nc, "longitude"); lat <- ncvar_get(nc, "latitude")
nc_close(nc)

grid_dt <- CJ(lon=lon, lat=lat)[, cell_id := .I]
n_cells <- nrow(grid_dt)

pts    <- vect(mri[, .(tmp_id, easting, northing)],
               geom=c("easting","northing"), crs="EPSG:27700")
coords <- as.data.table(crds(project(pts, "EPSG:4326")))
setnames(coords, c("lon_p","lat_p"))
mri <- cbind(mri, coords)

mri[, cell_id := mapply(function(x, y)
      grid_dt$cell_id[which.min((grid_dt$lon - x)^2 + (grid_dt$lat - y)^2)],
      lon_p, lat_p)]
unique_cells <- unique(mri$cell_id)

# Helper: read a NetCDF file and return a daily cell x day temperature matrix (C)
read_daily <- function(f) {
  nc <- nc_open(f)
  t2m <- ncvar_get(nc, "t2m") - 273.15
  tm  <- as.Date(as_datetime(ncvar_get(nc, "valid_time"), tz="UTC"))
  nc_close(nc)
  ud  <- unique(tm)
  out <- sapply(ud, function(d)
           rowMeans(matrix(t2m[,, tm == d], nrow=n_cells), na.rm=TRUE))
  list(mat = out[unique_cells, , drop=FALSE], dates = ud)
}

###############################################################
## STEP 4. Day-of-year baseline mean and cell-specific anomaly SD (1961-1990)
##   baseline mean       = mean of the 30 reference-period values per (cell, doy)
##   anomaly SD (sigma_i)= SD of the reference-period anomaly distribution per cell
###############################################################

files_baseline <- sort(list.files(nc_dir_baseline, pattern="\\.nc$", full.names=TRUE))
agg <- rbindlist(lapply(files_baseline, function(f) {
  d  <- read_daily(f)
  dt <- data.table(cell_id = rep(unique_cells, ncol(d$mat)),
                   date    = rep(d$dates, each=length(unique_cells)),
                   t2m     = as.vector(d$mat))
  dt[, md := format(date, "%m-%d")]
  dt <- dt[md != "02-29"][, doy := match(md, md_levels)]
  dt[, .(n=sum(!is.na(t2m)), s=sum(t2m, na.rm=TRUE), ss=sum(t2m^2, na.rm=TRUE)),
     by=.(cell_id, doy)]
}))[, .(n=sum(n), s=sum(s), ss=sum(ss)), by=.(cell_id, doy)]

baseline_dt <- agg[, .(cell_id, doy, mean_t2m_normal = s/n)]
# per-cell SD = sqrt( sum of squared deviations / (N - 1) )
sd_dcnt <- agg[, .(sd_dcnt = sqrt(sum(ss - s^2/n, na.rm=TRUE) /
                                  (sum(n, na.rm=TRUE) - 1))), by=cell_id]

# Coastal cells with no land baseline -> remap participants to nearest valid cell
na_cells <- baseline_dt[is.na(mean_t2m_normal), unique(cell_id)]
if (length(na_cells) > 0) {
  vc <- grid_dt[cell_id %in% baseline_dt[!is.na(mean_t2m_normal), unique(cell_id)]]
  remap <- sapply(na_cells, function(cid) {
    p <- grid_dt[cell_id == cid]
    vc$cell_id[which.min((vc$lon - p$lon)^2 + (vc$lat - p$lat)^2)]
  })
  mri[cell_id %in% na_cells, cell_id := remap[match(cell_id, na_cells)]]
}

###############################################################
## STEP 5. Daily anomaly and SD-based classification
###############################################################

files_obs <- sort(list.files(nc_dir_obs, pattern="\\.nc$", full.names=TRUE))
dcnt_daily <- rbindlist(lapply(files_obs, function(f) {
  yr <- as.integer(sub(".*t2m_(\\d{4})_\\d{2}\\.nc", "\\1", basename(f)))
  d  <- read_daily(f)
  dt <- data.table(cell_id = rep(unique_cells, ncol(d$mat)),
                   date    = rep(d$dates, each=length(unique_cells)),
                   t2m_obs = as.vector(d$mat), year = yr)
  dt[, md := format(date, "%m-%d")][md == "02-29", md := "02-28"]
  dt[, doy := match(md, md_levels)]
  dt <- merge(dt, baseline_dt, by=c("cell_id","doy"), all.x=TRUE)
  dt[, .(cell_id, year, dcnt = t2m_obs - mean_t2m_normal)]
}))

dcnt_daily <- merge(dcnt_daily, sd_dcnt, by="cell_id", all.x=TRUE)
dcnt_daily[, `:=`(
  warm_moderate = as.integer(dcnt >    sd_dcnt),
  warm_extreme  = as.integer(dcnt >  2*sd_dcnt),
  cold_moderate = as.integer(dcnt <   -sd_dcnt),
  cold_extreme  = as.integer(dcnt < -2*sd_dcnt))]

###############################################################
## STEP 6. Annual aggregation per cell
###############################################################

dcnt_annual <- dcnt_daily[, .(
  mean_dcnt     = mean(dcnt,         na.rm=TRUE),
  warm_mod_days = sum(warm_moderate, na.rm=TRUE),
  warm_ext_days = sum(warm_extreme,  na.rm=TRUE),
  cold_mod_days = sum(cold_moderate, na.rm=TRUE),
  cold_ext_days = sum(cold_extreme,  na.rm=TRUE)
), by=.(cell_id, year)]

###############################################################
## STEP 7. Cumulative anomaly-day exposure per participant (1-3 yr)
###############################################################

sum_vars <- c("warm_mod_days","warm_ext_days","cold_mod_days","cold_ext_days")
calc_exposure <- function(n_years) {
  m <- merge(mri[, .(tmp_id, cell_id, t0_year)], dcnt_annual,
             by="cell_id", allow.cartesian=TRUE)
  w <- m[year >= (t0_year - n_years) & year < t0_year]
  res <- w[, c(lapply(.SD[, ..sum_vars], sum, na.rm=TRUE),
               .(mean_dcnt = mean(mean_dcnt, na.rm=TRUE))), by=.(tmp_id)]
  setnames(res, c(sum_vars, "mean_dcnt"),
           paste0(c(sum_vars, "mean_dcnt"), "_", n_years, "yr"))
  res
}
exposure <- Reduce(function(a, b) merge(a, b, by="tmp_id"),
                   lapply(1:3, calc_exposure))

###############################################################
## STEP 8. Merge and save
###############################################################

final_dt <- merge(mri, exposure, by="tmp_id")
fwrite(final_dt, out_path)
cat("Saved:", out_path, "|", nrow(final_dt), "participants\n")
