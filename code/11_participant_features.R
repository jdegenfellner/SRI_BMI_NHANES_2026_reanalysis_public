# =====================================================================
# 11_participant_features.R  ->  participant_ndays.{rds,csv}
# Per-participant summaries of the minute-level series that the downstream
# figures need (number of valid days, missingness, and the day-pairing
# structure). Shipping this small table means the figures that characterise the
# disagreement (81, 83) reproduce WITHOUT the multi-hundred-MB minute series.
# Part of the raw pipeline: run after 10_build_common_series.R.
# =====================================================================
suppressMessages({library(data.table)})
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())
source(file.path(base, if (dir.exists(file.path(base, "code"))) "code" else "01_recompute_SRI", "00_paths.R"))

ser <- rbindlist(list(
  as.data.table(readRDS(file.path(inter, "series_1112.rds"))),
  as.data.table(readRDS(file.path(inter, "series_1314.rds")))))

# one row per participant: valid days, missingness, and the number of
# consecutive (adjacent day-of-week) vs non-consecutive day boundaries
feat <- ser[, {
  dw <- as.numeric(PAXDAYWM[seq(1, .N, by = 1440)])          # day-of-week per calendar day
  adj <- diff(dw) == 1 | (head(dw, -1) == 7 & tail(dw, -1) == 1)   # Sun->Mon wraps
  .(n_days = uniqueN(PAXDAYM),
    miss   = mean(is.na(binary_sleep)),
    n_adj  = sum(adj),
    n_gaps = sum(!adj))
}, by = SEQN]

saveRDS(feat, file.path(inter, "participant_ndays.rds"))
fwrite(feat, file.path(inter, "participant_ndays.csv"))
cat("wrote participant_ndays.{rds,csv}:", nrow(feat), "participants\n")
print(feat[, .(med_days = median(n_days), share_with_gap = round(mean(n_gaps > 0), 3))])
