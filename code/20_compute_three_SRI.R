# =====================================================================
# PART 2: Compute SRI three ways on the IDENTICAL NHANES sleep/wake
# series: (1) original hand formula, (2) sleepreg (real source fns),
# (3) GGIR::CalcSleepRegularityIndex.  All with package DEFAULT settings.
# =====================================================================
suppressMessages({library(data.table); library(GGIR)})
options(width = 160)

base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())
rdir   <- file.path(base, "01_recompute_SRI")
inter  <- file.path(rdir, "intermediate")

# ---- load real sleepreg functions (GGIR-free path) + percentile table
source(file.path(rdir, "sleepreg_src", "sleepreg_allFunctions.R"))
load(file.path(rdir, "sleepreg_src", "sysdata.rda"))   # -> quant (used by SRI_from_binary)

# analytic sample = the published main modeling frame (n=7085)
an_seqn <- unique(as.data.table(readRDS(file.path(base,
  "00_original_paper_materials/data/df_an1_imputed_MAIN.RDS")))$SEQN)
# load the two cycle series, restrict to analytic SEQN immediately (memory-light)
s1 <- as.data.table(readRDS(file.path(inter, "series_1112.rds")))[SEQN %in% an_seqn]
s2 <- as.data.table(readRDS(file.path(inter, "series_1314.rds")))[SEQN %in% an_seqn]
ser <- rbind(s1, s2); rm(s1, s2); gc()
setorder(ser, SEQN, UNIX_TIMESTAMP)
setkey(ser, SEQN)                       # binary-search subsetting in loops
seqns <- sort(unique(ser$SEQN))
cat("analytic SEQN:", length(an_seqn), " | present in filtered series:", length(seqns), "\n")

# ============================================================
# (1) HAND SRI  — replicate original calculate_sri (no imputation)
#     1-min epochs; weekday-consecutive day pairs (PAXDAYWM i vs i+1)
# ============================================================
hand_sri <- function(s, d, M = 1440L) {
  n <- length(s); if (n <= M) return(NA_real_)
  i   <- 1:(n - M)
  cur <- d[i]; nxt <- d[i + M]
  consec <- (nxt - cur == 1) | (cur == 7 & nxt == 1)
  a <- s[i]; b <- s[i + M]
  valid <- consec & !is.na(a) & !is.na(b)
  Nv <- sum(valid); if (Nv == 0) return(NA_real_)
  -100 + 200 / Nv * sum(a[valid] == b[valid])
}

# ============================================================
# (3) GGIR SRI — CalcSleepRegularityIndex per day-pair, mean over pairs
#     1-min epochs; consecutive CALENDAR dates, time-of-day matched
# ============================================================
ggir_sri <- function(ts, s) {
  df <- data.frame(
    time       = as.POSIXct(ts, origin = "1970-01-01", tz = "UTC"),
    sleepstate = as.numeric(s),
    invalid    = as.integer(is.na(s))
  )
  out <- tryCatch(
    GGIR::CalcSleepRegularityIndex(data = df, epochsize = 60, desiredtz = "UTC"),
    error = function(e) NA)
  if (is.data.frame(out)) mean(out$SleepRegularityIndex, na.rm = TRUE) else NA_real_
}

# ---- compute hand + GGIR in one pass over participants ----
res <- data.table(SEQN = seqns, SRI_hand = NA_real_, SRI_GGIR = NA_real_)
t0 <- Sys.time()
for (k in seq_along(seqns)) {
  d <- ser[.(seqns[k])]
  res$SRI_hand[k] <- hand_sri(d$binary_sleep, d$PAXDAYWM)
  res$SRI_GGIR[k] <- ggir_sri(d$UNIX_TIMESTAMP, d$binary_sleep)
  if (k %% 500 == 0) cat(sprintf("  hand/GGIR %d/%d  (%.1f min)\n",
        k, length(seqns), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
saveRDS(res, file.path(inter, "sri_hand_ggir.rds"))
cat("hand+GGIR done\n")

# ============================================================
# (2) sleepreg SRI — real SRI_from_binary (10-s epochs, continuous
#     24h-lag, defaults minSRIdays=5, exclNAhrs=6).  File-based.
# ============================================================
bindir <- file.path(inter, "sleepreg_binary_in")
unlink(bindir, recursive = TRUE); dir.create(bindir, showWarnings = FALSE)
# write 2-col CSV per participant: col1 binary_sleep, col2 UNIX timestamp
for (q in seqns) {
  d <- ser[.(q)]
  fwrite(data.table(binary_sleep = d$binary_sleep, ts = d$UNIX_TIMESTAMP),
         file.path(bindir, paste0("sw_", q, ".csv")))
}
cat("wrote", length(seqns), "binary input files\n")
# minute series -> transition (SWS) format
SWS_from_binarySW(binarySWdir = bindir, tsCol = 2, binaryCol = 1)
swsdir <- file.path(bindir, "SWS_output")
# SRI from SWS files (package defaults; no raster for speed)
alloutdir <- file.path(inter, "sleepreg_out")
unlink(alloutdir, recursive = TRUE); dir.create(alloutdir, showWarnings = FALSE)
SRI_from_binary(binarydir = swsdir, alloutdir = alloutdir,
                col.trans = 1, col.timestamp = 2,
                wr.raster = FALSE, overwr = TRUE)   # defaults: minSRIdays=5, exclNAhrs=6
sr <- fread(file.path(alloutdir, "SRI.csv"))
sr[, SEQN := as.numeric(gsub("\\D", "", File))]
sr <- sr[, .(SEQN, SRI_sleepreg = as.numeric(SRI), sleepreg_days = as.numeric(SRI_days))]
saveRDS(sr, file.path(inter, "sri_sleepreg.rds"))
cat("sleepreg done\n")

# ============================================================
# Merge all three + original published SRI; quick comparison
# ============================================================
res <- readRDS(file.path(inter, "sri_hand_ggir.rds"))
all <- merge(res, sr, by = "SEQN", all = TRUE)
an  <- as.data.table(readRDS(file.path(base,
        "00_original_paper_materials/data/Combined_data_sets_2_waves_analytic_data_set.RDS")))
an  <- an[, .(SEQN, SRI_orig = SRI_no_imputation, BMXBMI, RIAGENDR, RIDRETH3,
              RIDAGEYR, WTMEC2YR, SDMVPSU, SDMVSTRA, SDDSRVYR)]
all <- merge(all, an, by = "SEQN", all.x = TRUE)
saveRDS(all, file.path(inter, "sri_all_methods.rds"))
fwrite(all, file.path(inter, "sri_all_methods.csv"))

cat("\n=== VALIDATION: hand vs original SRI_no_imputation (analytic only) ===\n")
v <- all[!is.na(SRI_orig) & !is.na(SRI_hand)]
cat("n compared:", nrow(v), "\n")
cat("Pearson r:", round(cor(v$SRI_hand, v$SRI_orig), 5), "\n")
cat("max |diff|:", round(max(abs(v$SRI_hand - v$SRI_orig)), 4),
    " mean |diff|:", round(mean(abs(v$SRI_hand - v$SRI_orig)), 4), "\n")

cat("\n=== SRI summaries (analytic sample) ===\n")
a <- all[!is.na(SRI_orig)]
for (col in c("SRI_orig","SRI_hand","SRI_GGIR","SRI_sleepreg")) {
  x <- a[[col]]
  cat(sprintf("%-14s n=%4d  mean=%6.2f sd=%5.2f  median=%6.2f  [%.1f, %.1f]\n",
      col, sum(!is.na(x)), mean(x,na.rm=T), sd(x,na.rm=T), median(x,na.rm=T),
      min(x,na.rm=T), max(x,na.rm=T)))
}
cat("\n=== pairwise correlations (complete cases, analytic) ===\n")
print(round(cor(a[, .(SRI_orig, SRI_hand, SRI_GGIR, SRI_sleepreg)], use = "pairwise.complete.obs"), 3))
cat("\nDONE part 2\n")
