###############################################################
## 01_data_preparation.R
##
## Build the participant-level analysis dataset.
##  - ERA5-Land 1961-1990 day-of-year baseline mean (mu_id) and
##    cell-specific SD of daily anomalies (sigma_i)
##  - daily temperature anomaly = T - mu_id; classification into
##    warm/cold (>±1 SD) and extreme warm/cold (>±2 SD) anomaly days
##  - annual and cumulative (1-3 yr) anomaly-day counts per participant
##  - longitudinal brain volume change (t1 - t0) and covariates
##
## Output: data/mri_dcnt_daily.csv
###############################################################
library(terra)
library(ncdf4)
library(data.table)
library(lubridate)

###############################################################
## 경로 (Windows)
###############################################################

BASE        <- "."
mri_path    <- file.path(BASE, "3.R/data/data_brain2.csv")
era5_sample <- file.path(BASE, "3.R/data/heatwave_final/t2m_2014_01.nc")
cell_annual_path <- file.path(BASE, "3.R/data/heatwave_final/exposure_by_cell_year.csv")
nc_dir_baseline  <- file.path(BASE, "data/ERA5/unzipped_1961_1990")
nc_dir_obs       <- file.path(BASE, "data/ERA5/unzipped_real")
dcnt_dir         <- file.path(BASE, "data/ERA5/dcnt")
out_path         <- file.path(BASE, "3.R/data/mri_dcnt_daily.csv")   # ★ _daily

dir.create(dcnt_dir, showWarnings=FALSE)

## ★ 파라미터
MEAN_WINDOW <- 0        # 0 = 해당 캘린더일 30개만(창 없음); >0 이면 +-N일 이동창 평균
SD_PERIOD   <- "base"   # "base" = 1961~1990 기준기간 dcnt 분포 SD (기본); "obs" = 관측기간 dcnt
md_levels   <- format(seq(as.Date("2001-01-01"), as.Date("2001-12-31"), by="day"), "%m-%d")

###############################################################
## STEP 1. MRI 로드 & 뇌 영상 변화량 계산
###############################################################

cat("===== STEP 1. MRI 로드 & 변화량 계산 =====\n")

mri <- fread(mri_path)
cat("원자료:", nrow(mri), "명\n")

mri[, t0_date := as.Date(p53_i2)]
mri[, t1_date := as.Date(p53_i3)]
mri[, tmp_id  := .I]

mri_valid <- mri[
  !is.na(t0_date) & !is.na(t1_date) &
    !is.na(p25009_i2) & !is.na(p25009_i3) &
    !is.na(p32223_a0) & !is.na(p32224_a0)
]

mri_valid[, t0_year       := year(t0_date)]
mri_valid[, t1_year       := year(t1_date)]
mri_valid[, followup_years := as.numeric(t1_date - t0_date) / 365.25]
cat("분석 가능:", nrow(mri_valid), "명\n")

mri_valid[, p25781_i2_log := log(as.numeric(p25781_i2) + 1)]
mri_valid[, p25781_i3_log := log(as.numeric(p25781_i3) + 1)]

brain_pairs <- list(
  c("p25009_i2","p25009_i3","delta_brain_norm"),
  c("p25005_i2","p25005_i3","delta_gm_norm"),
  c("p25007_i2","p25007_i3","delta_wm_norm"),
  c("p25010_i2","p25010_i3","delta_brain"),
  c("p25006_i2","p25006_i3","delta_gm"),
  c("p25008_i2","p25008_i3","delta_wm"),
  c("p25781_i2_log","p25781_i3_log","delta_wmh_log"),
  c("p25019_i2","p25019_i3","delta_hipp_L"),
  c("p25020_i2","p25020_i3","delta_hipp_R")
)
for (p in brain_pairs)
  mri_valid[, (p[3]) := as.numeric(get(p[2])) - as.numeric(get(p[1]))]

mri_valid[, delta_hipp_bilateral :=
            (as.numeric(p25019_i3) + as.numeric(p25020_i3)) -
            (as.numeric(p25019_i2) + as.numeric(p25020_i2))]

delta_cols <- grep("^delta_", names(mri_valid), value=TRUE)
for (d in delta_cols)
  mri_valid[, (sub("^delta_","rate_",d)) := get(d) / followup_years]

cat("변화량 계산 완료\n")

###############################################################
## STEP 2. t0 기준 주소 찾기 (이사 이력 반영)
###############################################################

cat("\n===== STEP 2. 주소 이력 및 좌표 추출 =====\n")

mri_valid[, best_array := 0L]
mri_valid[!is.na(p32220_a1) & as.Date(p32220_a1) <= t0_date, best_array := 1L]
mri_valid[!is.na(p32220_a2) & as.Date(p32220_a2) <= t0_date, best_array := 2L]
mri_valid[!is.na(p32220_a3) & as.Date(p32220_a3) <= t0_date, best_array := 3L]
mri_valid[!is.na(p32220_a4) & as.Date(p32220_a4) <= t0_date, best_array := 4L]

cat("이사 이력 분포:\n"); print(table(mri_valid$best_array))

mri_valid[, easting := as.integer(mapply(
  function(arr, ...) c(...)[[arr+1]],
  best_array, p32223_a0, p32223_a1, p32223_a2, p32223_a3, p32223_a4))]
mri_valid[, northing := as.integer(mapply(
  function(arr, ...) c(...)[[arr+1]],
  best_array, p32224_a0, p32224_a1, p32224_a2, p32224_a3, p32224_a4))]

mri_valid <- mri_valid[!is.na(easting) & !is.na(northing)]
cat("좌표 있는 분석 가능:", nrow(mri_valid), "명\n")

###############################################################
## STEP 3. ERA5 grid cell 매칭 + 연간 노출 계산
###############################################################

cat("\n===== STEP 3. ERA5 cell 매칭 + 연간 노출 =====\n")

nc  <- nc_open(era5_sample)
lon <- ncvar_get(nc, "longitude")
lat <- ncvar_get(nc, "latitude")
nc_close(nc)

grid_dt <- CJ(lon=lon, lat=lat)
grid_dt[, cell_id := .I]
n_cells <- nrow(grid_dt)

pts     <- vect(mri_valid[, .(tmp_id, easting, northing)],
                geom=c("easting","northing"), crs="EPSG:27700")
pts_wgs <- project(pts, "EPSG:4326")
coords  <- as.data.table(crds(pts_wgs))
setnames(coords, c("lon_p","lat_p"))
mri_valid <- cbind(mri_valid, coords)

find_cell <- function(lon_p, lat_p) {
  dist <- (grid_dt$lon - lon_p)^2 + (grid_dt$lat - lat_p)^2
  grid_dt$cell_id[which.min(dist)]
}
cat("cell_id 매칭 중...\n")
mri_valid[, cell_id := mapply(find_cell, lon_p, lat_p)]
cat("고유 cell 수:", uniqueN(mri_valid$cell_id), "\n")

unique_cells <- unique(mri_valid$cell_id)

cell_annual <- fread(cell_annual_path)

dt_all <- merge(
  mri_valid[, .(tmp_id, cell_id, t0_year, t1_year, followup_years)],
  cell_annual, by="cell_id", allow.cartesian=TRUE
)

exposure_vars_sum  <- c("heat_days_95","heat_days_97.5","heat_days_99",
                        "cold_days_5","cold_days_2.5","cold_days_1")
exposure_vars_mean <- c("mean_t2m","mean_rh")

calc_exposure <- function(dt, label, year_filter_expr) {
  cat("  계산 창:", label, "\n")
  dt_w   <- dt[eval(year_filter_expr)]
  result <- dt_w[, c(
    lapply(.SD[, ..exposure_vars_sum],  sum,  na.rm=TRUE),
    lapply(.SD[, ..exposure_vars_mean], mean, na.rm=TRUE)
  ), by=.(tmp_id, t0_year, t1_year, followup_years)]
  setnames(result,
           c(exposure_vars_sum, exposure_vars_mean),
           paste0(c(exposure_vars_sum, exposure_vars_mean), "_", label))
  result
}

exp_windows <- list(
  calc_exposure(dt_all,"1yr",        quote(year >= (t0_year-1) & year < t0_year)),
  calc_exposure(dt_all,"2yr",        quote(year >= (t0_year-2) & year < t0_year)),
  calc_exposure(dt_all,"3yr",        quote(year >= (t0_year-3) & year < t0_year)),
  calc_exposure(dt_all,"4yr",        quote(year >= (t0_year-4) & year < t0_year)),
  calc_exposure(dt_all,"5yr",        quote(year >= (t0_year-5) & year < t0_year)),
  calc_exposure(dt_all,"cumulative", quote(year >= t0_year & year <= t1_year))
)

exposure_dt <- Reduce(
  function(a,b) merge(a, b, by=c("tmp_id","t0_year","t1_year","followup_years")),
  exp_windows
)
cat("연간 노출 계산 완료:", nrow(exposure_dt), "명\n")

mri_valid[, c("lon_p","lat_p") := NULL]
final_dt <- merge(mri_valid, exposure_dt,
                  by=c("tmp_id","t0_year","t1_year","followup_years"))
cat("연간 노출 merge 완료:", nrow(final_dt), "명\n")

###############################################################
## STEP 4. doy 기준선 평균 + 셀별 dcnt SD 산출 (1961~1990)   ★ 변경
##   - 각 (cell, doy) 의 기준선 평균 = 해당 셀, 해당 캘린더일 30개 평균
##   - 셀별 dcnt SD = 기준기간(1961~1990) 이상치 분포의 표준편차
##     (각 doy 의 제곱편차합 ss - s^2/n 을 셀 단위로 합산)
##   - MEAN_WINDOW > 0 이면 평균은 +-N일 이동창 (단 이때 SD 는 doy 평균 기준)
###############################################################

cat("\n===== STEP 4. doy 기준선 평균 + 기준기간 dcnt SD =====\n")

baseline_path <- file.path(dcnt_dir, "baseline_1961_1990_daily.csv")
sd_dcnt_path  <- file.path(dcnt_dir, "sd_dcnt_daily.csv")

if (file.exists(baseline_path) && file.exists(sd_dcnt_path)) {
  cat("기존 daily 기준선 평균 + sd_dcnt 로드\n")
  baseline_dt <- fread(baseline_path)
  sd_dcnt     <- fread(sd_dcnt_path)

} else {
  files_baseline <- sort(list.files(nc_dir_baseline, pattern="\\.nc$", full.names=TRUE))
  cat("baseline 파일 수:", length(files_baseline), "개\n")

  agg_list <- vector("list", length(files_baseline))
  for (i in seq_along(files_baseline)) {
    nc       <- nc_open(files_baseline[i])
    t2m_raw  <- ncvar_get(nc, "t2m")
    time_raw <- ncvar_get(nc, "valid_time")
    nc_close(nc)

    t2m_c <- t2m_raw - 273.15
    dates <- as.Date(as_datetime(time_raw, tz="UTC"))
    ud    <- unique(dates)

    daily_t2m <- matrix(NA, nrow=n_cells, ncol=length(ud))
    for (d in seq_along(ud)) {
      idx <- which(dates == ud[d])
      daily_t2m[,d] <- rowMeans(matrix(t2m_c[,,idx], nrow=n_cells), na.rm=TRUE)
    }

    ct  <- daily_t2m[unique_cells,, drop=FALSE]
    dtf <- data.table(
      cell_id = rep(unique_cells, ncol(ct)),
      date    = rep(ud, each=length(unique_cells)),
      t2m     = as.vector(ct)
    )
    dtf[, md := format(date, "%m-%d")]
    dtf <- dtf[md != "02-29"]
    dtf[, doy := match(md, md_levels)]
    agg_list[[i]] <- dtf[, .(n=sum(!is.na(t2m)), s=sum(t2m, na.rm=TRUE),
                             ss=sum(t2m^2, na.rm=TRUE)),
                         by=.(cell_id, doy)]
    rm(t2m_raw, t2m_c, daily_t2m, ct, dtf); gc()
  }
  agg <- rbindlist(agg_list)[, .(n=sum(n), s=sum(s), ss=sum(ss)),
                             by=.(cell_id, doy)]

  if (MEAN_WINDOW > 0) {
    # +-MEAN_WINDOW 순환 이동창 평균
    base_mean <- NULL
    for (o in -MEAN_WINDOW:MEAN_WINDOW) {
      tmp <- agg[, .(cell_id, doy=((doy-1+o) %% 365)+1, n, s)]
      base_mean <- if (is.null(base_mean)) tmp else
        rbind(base_mean, tmp)[, .(n=sum(n), s=sum(s)), by=.(cell_id, doy)]
    }
    baseline_dt <- base_mean[, .(cell_id, doy, mean_t2m_normal = s/n)]
  } else {
    baseline_dt <- agg[, .(cell_id, doy, mean_t2m_normal = s/n)]
  }
  fwrite(baseline_dt, baseline_path)

  # ★ 셀별 기준기간 dcnt 분포 SD
  #   각 doy 의 제곱편차합 = ss - s^2/n; 셀 단위 합산 후 sqrt(SSD/(Ntot-1))
  agg[, ssd := ss - s^2 / n]
  sd_dcnt <- agg[, .(sd_dcnt = sqrt(sum(ssd, na.rm=TRUE) /
                                    (sum(n, na.rm=TRUE) - 1))), by=cell_id]
  fwrite(sd_dcnt, sd_dcnt_path)
}
cat("기준선 평균 완료:", nrow(baseline_dt), "행 |",
    "cell:", uniqueN(baseline_dt$cell_id), "|",
    "평균 sd_dcnt:", round(mean(sd_dcnt$sd_dcnt, na.rm=TRUE), 2), "C\n")

###############################################################
## STEP 4b. 해안(NA) cell 사전 보정
##   - 바다/해안 격자(기준선 NaN)에 매핑된 참가자를
##     가장 가까운 유효(육지) cell 로 노출 계산 전에 재매핑
##   - 따라서 exp_dcnt_all 부터 0/NaN 이 생기지 않음
###############################################################

cat("\n===== STEP 4b. 해안 cell 사전 보정 =====\n")

na_cells_baseline <- baseline_dt[is.na(mean_t2m_normal), unique(cell_id)]
cat("기준선 NA(바다) cell 수:", length(na_cells_baseline), "\n")

if (length(na_cells_baseline) > 0) {
  valid_cells  <- baseline_dt[!is.na(mean_t2m_normal), unique(cell_id)]
  valid_coords <- grid_dt[cell_id %in% valid_cells]

  find_nearest_valid <- function(cid) {
    na_coord <- grid_dt[cell_id == cid, .(lon, lat)]
    d <- (valid_coords$lon - na_coord$lon)^2 + (valid_coords$lat - na_coord$lat)^2
    valid_coords$cell_id[which.min(d)]
  }

  remap <- data.table(
    cell_id_orig    = na_cells_baseline,
    cell_id_replace = sapply(na_cells_baseline, find_nearest_valid)
  )

  n_fixed <- final_dt[cell_id %in% na_cells_baseline, .N]
  final_dt[cell_id %in% na_cells_baseline,
           cell_id := remap$cell_id_replace[match(cell_id, remap$cell_id_orig)]]
  cat("재매핑된 참가자:", n_fixed, "명 (가장 가까운 유효 cell 로)\n")
}

###############################################################
## STEP 5. DCNT 계산 + 셀별 dcnt SD 로 극한 분류        ★ 변경
###############################################################

cat("\n===== STEP 5. DCNT 계산 + dcnt SD 분류 =====\n")

dcnt_daily_path <- file.path(dcnt_dir, "dcnt_daily_daily.csv")

if (file.exists(dcnt_daily_path)) {
  cat("기존 daily DCNT(원시) 로드\n")
  dcnt_daily <- fread(dcnt_daily_path)

} else {
  files_obs <- sort(list.files(nc_dir_obs, pattern="\\.nc$", full.names=TRUE))
  cat("관측 파일 수:", length(files_obs), "개\n")

  # 5-1. 모든 관측일의 dcnt 계산 (분류는 SD 산출 후)
  dcnt_list <- vector("list", length(files_obs))
  for (i in seq_along(files_obs)) {
    f  <- files_obs[i]
    yr <- as.integer(sub(".*t2m_(\\d{4})_\\d{2}\\.nc","\\1", basename(f)))
    mo <- as.integer(sub(".*t2m_\\d{4}_(\\d{2})\\.nc","\\1", basename(f)))
    cat("  처리:", yr, "/", sprintf("%02d", mo), "\n")

    nc       <- nc_open(f)
    t2m_raw  <- ncvar_get(nc, "t2m")
    time_raw <- ncvar_get(nc, "valid_time")
    nc_close(nc)

    t2m_c <- t2m_raw - 273.15
    dates <- as.Date(as_datetime(time_raw, tz="UTC"))
    ud    <- unique(dates)

    daily_t2m <- matrix(NA, nrow=n_cells, ncol=length(ud))
    for (d in seq_along(ud)) {
      idx <- which(dates == ud[d])
      daily_t2m[,d] <- rowMeans(matrix(t2m_c[,,idx], nrow=n_cells), na.rm=TRUE)
    }

    ct <- daily_t2m[unique_cells,, drop=FALSE]
    dt_nc <- data.table(
      cell_id = rep(unique_cells, ncol(ct)),
      date    = rep(ud, each=length(unique_cells)),
      t2m_obs = as.vector(ct),
      year=yr, month=mo
    )
    dt_nc[, md := format(date, "%m-%d")]
    dt_nc[md == "02-29", md := "02-28"]
    dt_nc[, doy := match(md, md_levels)]
    dt_nc <- merge(dt_nc, baseline_dt, by=c("cell_id","doy"), all.x=TRUE)
    dt_nc[, dcnt := t2m_obs - mean_t2m_normal]

    dcnt_list[[i]] <- dt_nc[, .(cell_id, date, year, month, doy, t2m_obs, dcnt)]
    rm(t2m_raw, t2m_c, daily_t2m, ct, dt_nc); gc()
  }
  dcnt_daily <- rbindlist(dcnt_list)
  fwrite(dcnt_daily, dcnt_daily_path)   # 원시 dcnt 저장 (분류 플래그는 아래에서 in-memory)
}

# 5-2. SD_PERIOD=="obs" 면 관측기간 dcnt 분포로 SD 재산출 (기본은 STEP 4 의 기준기간 SD)
if (SD_PERIOD == "obs") {
  sd_dcnt <- dcnt_daily[, .(sd_dcnt = sd(dcnt, na.rm=TRUE)), by=cell_id]
  fwrite(sd_dcnt, sd_dcnt_path)
}

# 5-3. 셀별 dcnt SD 로 극한 분류
if ("sd_dcnt" %in% names(dcnt_daily)) dcnt_daily[, sd_dcnt := NULL]
dcnt_daily <- merge(dcnt_daily, sd_dcnt, by="cell_id", all.x=TRUE)
dcnt_daily[, warm_moderate := as.integer(dcnt >  sd_dcnt)]
dcnt_daily[, warm_extreme  := as.integer(dcnt >  2*sd_dcnt)]
dcnt_daily[, cold_moderate := as.integer(dcnt < -sd_dcnt)]
dcnt_daily[, cold_extreme  := as.integer(dcnt < -2*sd_dcnt)]

cat("DCNT 계산 완료:", nrow(dcnt_daily), "행 | SD_PERIOD =", SD_PERIOD, "|",
    "평균 sd_dcnt:", round(mean(sd_dcnt$sd_dcnt, na.rm=TRUE), 2), "C\n")

###############################################################
## STEP 6. DCNT 연간 집계
###############################################################

cat("\n===== STEP 6. 연간 집계 =====\n")

dcnt_annual_path <- file.path(dcnt_dir, "dcnt_annual_daily.csv")

if (file.exists(dcnt_annual_path)) {
  cat("기존 daily 연간 집계 로드\n")
  dcnt_annual <- fread(dcnt_annual_path)
} else {
  dcnt_annual <- dcnt_daily[, .(
    mean_dcnt     = mean(dcnt,          na.rm=TRUE),
    warm_mod_days = sum(warm_moderate,  na.rm=TRUE),
    warm_ext_days = sum(warm_extreme,   na.rm=TRUE),
    cold_mod_days = sum(cold_moderate,  na.rm=TRUE),
    cold_ext_days = sum(cold_extreme,   na.rm=TRUE)
  ), by=.(cell_id, year)]
  fwrite(dcnt_annual, dcnt_annual_path)
}
cat("연간 집계 완료:", nrow(dcnt_annual), "행\n")

###############################################################
## STEP 7. 참가자별 DCNT 노출 계산 (1~5yr + cumulative)
###############################################################

cat("\n===== STEP 7. 참가자별 DCNT 노출 계산 =====\n")

dcnt_vars_sum  <- c("warm_mod_days","warm_ext_days","cold_mod_days","cold_ext_days")
dcnt_vars_mean <- c("mean_dcnt")

calc_dcnt_exposure <- function(n_years, label) {
  cat("  계산 창:", label, "\n")
  dt_m <- merge(
    final_dt[, .(tmp_id, cell_id, t0_year)],
    dcnt_annual, by="cell_id", allow.cartesian=TRUE
  )
  dt_w <- dt_m[year >= (t0_year - n_years) & year < t0_year]
  res  <- dt_w[, c(
    lapply(.SD[, ..dcnt_vars_sum],  sum,  na.rm=TRUE),
    lapply(.SD[, ..dcnt_vars_mean], mean, na.rm=TRUE)
  ), by=.(tmp_id)]
  setnames(res,
           c(dcnt_vars_sum, dcnt_vars_mean),
           paste0(c(dcnt_vars_sum, dcnt_vars_mean), "_", label))
  res
}

exp_list <- lapply(1:5, function(n) calc_dcnt_exposure(n, paste0(n,"yr")))

cat("  계산 창: cumulative\n")
dt_cum <- merge(
  final_dt[, .(tmp_id, cell_id, t0_year, t1_year)],
  dcnt_annual, by="cell_id", allow.cartesian=TRUE
)
exp_cum <- dt_cum[year >= t0_year & year <= t1_year][, c(
  lapply(.SD[, ..dcnt_vars_sum],  sum,  na.rm=TRUE),
  lapply(.SD[, ..dcnt_vars_mean], mean, na.rm=TRUE)
), by=.(tmp_id)]
setnames(exp_cum,
         c(dcnt_vars_sum, dcnt_vars_mean),
         paste0(c(dcnt_vars_sum, dcnt_vars_mean), "_cumulative"))

exp_list[[6]] <- exp_cum
exp_dcnt_all  <- Reduce(function(a,b) merge(a,b, by="tmp_id"), exp_list)
cat("DCNT 노출 계산 완료\n")

###############################################################
## STEP 8. 전처리 + 전체 merge
###############################################################

cat("\n===== STEP 8. 전처리 & merge =====\n")

final_dt_dcnt <- merge(final_dt, exp_dcnt_all, by="tmp_id")

na_like <- c("Prefer not to answer","Do not know","None of the above","")

final_dt_dcnt[, p31 := factor(p31, levels=c("Male","Female"))]

final_dt_dcnt[p20116_i2 %in% na_like, p20116_i2 := NA_character_]
final_dt_dcnt[, p20116_i2 := factor(p20116_i2,
                                    levels=c("Never","Previous","Current"))]

final_dt_dcnt[p1558_i2 %in% na_like, p1558_i2 := NA_character_]
final_dt_dcnt[, p1558_i2 := factor(p1558_i2,
                                   levels=c("Never","Special occasions only","One to three times a month",
                                            "Once or twice a week","Three or four times a week",
                                            "Daily or almost daily"))]

final_dt_dcnt[p2443_i2 %in% na_like, p2443_i2 := NA_character_]
final_dt_dcnt[, p2443_i2 := factor(
  ifelse(p2443_i2=="Yes","Yes","No"), levels=c("No","Yes"))]

raw_pa  <- final_dt_dcnt$p6164_i2
raw_pa[raw_pa %in% na_like] <- NA_character_
pipe_pa <- ifelse(is.na(raw_pa), NA_character_, paste0("|",raw_pa,"|"))
pa_choices <- c(
  "Walking for pleasure (not as a means of transport)",
  "Other exercises (eg: swimming, cycling, keep fit, bowling)",
  "Strenuous sports","Light DIY (eg: pruning, watering the lawn)",
  "Heavy DIY (eg: weeding, lawn mowing, carpentry, digging)")
pa_hits <- rowSums(sapply(pa_choices, function(ch)
  !is.na(pipe_pa) & grepl(paste0("|",ch,"|"), pipe_pa, fixed=TRUE)))
final_dt_dcnt[, pa_any := as.integer(pa_hits > 0)]

final_dt_dcnt[p21000_i0 %in% c("Prefer not to answer","Do not know"),
              p21000_i0 := NA_character_]
final_dt_dcnt[, ethnicity := fcase(
  p21000_i0 %in% c("White","Any other white background","British","Irish"), "White",
  p21000_i0 %in% c("Asian or Asian British","Any other Asian background",
                   "Indian","Pakistani","Bangladeshi","Chinese"),           "Asian",
  p21000_i0 %in% c("Black or Black British","Any other Black background",
                   "African","Caribbean"),                                  "Black",
  p21000_i0 %in% c("Mixed","Any other mixed background","White and Asian",
                   "White and Black African","White and Black Caribbean"),  "Mixed",
  p21000_i0 == "Other ethnic group",                                        "Other",
  default=NA_character_)]
final_dt_dcnt[, ethnicity := factor(ethnicity,
                                    levels=c("White","Asian","Black","Mixed","Other"))]

final_dt_dcnt[p6138_i0 %in% c("Prefer not to answer","Do not know"),
              p6138_i0 := NA_character_]
final_dt_dcnt[, edu_grouped := fcase(
  grepl("College or University degree",     p6138_i0, fixed=TRUE), "Degree",
  grepl("A levels/AS levels or equivalent", p6138_i0, fixed=TRUE) |
    grepl("Other professional qualifications",p6138_i0, fixed=TRUE), "A-level",
  grepl("O levels/GCSEs or equivalent",    p6138_i0, fixed=TRUE) |
    grepl("CSEs or equivalent",             p6138_i0, fixed=TRUE) |
    grepl("NVQ or HND or HNC or equivalent",p6138_i0, fixed=TRUE), "O-level",
  default=NA_character_)]
final_dt_dcnt[, edu_grouped := factor(edu_grouped,
                                      levels=c("O-level","A-level","Degree"))]

final_dt_dcnt[, p54_i2 := factor(p54_i2)]
final_dt_dcnt[, p54_i3 := factor(p54_i3)]

cat("\n표본:", nrow(final_dt_dcnt), "명 x", ncol(final_dt_dcnt), "변수\n")

###############################################################
## STEP 9. NA 점검 (보정은 STEP 4b 에서 사전 처리됨)
###############################################################

cat("\n===== STEP 9. NA 점검 =====\n")
cat("mean_dcnt_1yr NaN/NA:",   sum(is.na(final_dt_dcnt$mean_dcnt_1yr)),      "\n")
cat("warm_ext_days_1yr NA:",   sum(is.na(final_dt_dcnt$warm_ext_days_1yr)), "\n")
cat("warm_ext_days_1yr == 0:", sum(final_dt_dcnt$warm_ext_days_1yr == 0, na.rm=TRUE), "명 (정상 0 포함)\n")


###############################################################
## STEP 10. 최종 저장
###############################################################

cat("\n===== STEP 10. 최종 저장 =====\n")

fwrite(final_dt_dcnt, out_path)
cat("저장 완료:", out_path, "\n")
cat("\n부산물:\n")
cat("  -", baseline_path,    "(cell, doy, mean_t2m_normal)\n")
cat("  -", sd_dcnt_path,     "(cell, sd_dcnt)\n")
cat("  -", dcnt_daily_path,  "\n")
cat("  -", dcnt_annual_path, "\n")
