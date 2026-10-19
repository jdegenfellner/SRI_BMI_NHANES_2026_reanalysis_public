# =====================================================================
# Follow-up commentary: recompute SRI with sleepreg & GGIR vs hand calc
# PART 1: Rebuild the IDENTICAL filtered per-minute sleep/wake series
#         that the original paper used, with correct clock-anchored
#         UNIX timestamps. Output: one compact RDS per cycle + combined.
# =====================================================================
suppressMessages({library(data.table); library(haven)})
options(width = 160)

base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())
ddir   <- file.path(base, "01_recompute_SRI", "data_local")
outdir <- file.path(base, "01_recompute_SRI", "intermediate")
dir.create(outdir, showWarnings = FALSE)

# ---- helper: replicate original inclusion filters (steps 1-7) -------
filter_cycle <- function(pax, header_xpt, cycle_label) {
  setDT(pax)
  # keep only needed columns to save memory
  keep <- intersect(c("SEQN","PAXDAYM","PAXDAYWM","PAXSSNMP","PAXTSM","PAXPREDM","PAXQFM"), names(pax))
  pax  <- pax[, ..keep]
  pax[, PAXPREDM := as.character(PAXPREDM)]

  # Step 0: binary_sleep (sleep=1, wake=0, non-wear/unknown=NA)
  pax[, binary_sleep := fifelse(PAXPREDM=="sleep", 1, fifelse(PAXPREDM=="wake", 0, NA_real_))]
  cat(cycle_label, "initial SEQN:", uniqueN(pax$SEQN), "\n")

  # Step 1: full-minute observations
  pax <- pax[PAXTSM == 60]
  # Step 2: quality flag
  pax <- pax[PAXQFM == 0]
  cat("  after PAXTSM==60 & PAXQFM==0:", uniqueN(pax$SEQN), "SEQN\n")

  # Step 3: >=3 days with <120 min non-wear
  nw <- pax[, .(non_wear = sum(PAXPREDM=="non-wear")), by=.(SEQN, PAXDAYM)][non_wear < 120]
  good_nw <- nw[, .(vd = uniqueN(PAXDAYM)), by=SEQN][vd >= 3, SEQN]
  pax <- pax[SEQN %in% good_nw]
  cat("  after non-wear filter:", uniqueN(pax$SEQN), "SEQN\n")

  # Step 4: keep only days with exactly 1440 timesteps
  pax[, nts := .N, by=.(SEQN, PAXDAYM)]
  pax <- pax[nts == 1440][, nts := NULL]
  cat("  after 1440-timestep filter:", uniqueN(pax$SEQN), "SEQN\n")

  # Step 5: >=3 consecutive days (rle diff>=2) incl. weekend (PAXDAYWM 1 or 7)
  keep_seqn <- pax[, {
    ds <- sort(unique(as.numeric(PAXDAYWM)))
    consec <- if (sum(!is.na(ds)) < 3) FALSE else max(rle(diff(ds))$lengths) >= 2
    wknd   <- any(ds %in% c(1,7))
    .(ok = consec && wknd)
  }, by=SEQN][ok==TRUE, SEQN]
  pax <- pax[SEQN %in% keep_seqn]
  cat("  after consecutive+weekend filter:", uniqueN(pax$SEQN), "SEQN\n")

  # Step 6: missingness of binary_sleep < 30%
  miss <- pax[, .(mr = mean(is.na(binary_sleep))), by=SEQN][mr < 0.30, SEQN]
  pax <- pax[SEQN %in% miss]
  cat("  after <30% missingness filter:", uniqueN(pax$SEQN), "SEQN  (rows:", nrow(pax), ")\n")

  # ---- clock-anchored UNIX timestamp ----
  # Every retained day is a complete calendar day (exactly 1440 minutes, step 4), so its
  # first minute is 00:00:00 local time. Anchor each day on its own first sample; no header
  # start time is needed. Days are placed (PAXDAYM - 1) * 86400 s apart, so a dropped day
  # leaves a real gap in the timeline (day 1 = day of monitor placement).
  # Note (fixed 2026-10-14): the previous version derived the offset from PAXHD$PAXFTIME.
  # That column is an empty string for ~12% of participants, and as.POSIXct() on the whole
  # vector without an explicit format then silently falls back to a date-only format, so
  # every start time became 00:00:00 and every series was shifted by the placement time.
  pax[, UNIX_TIMESTAMP := (as.integer(PAXDAYM) - 1) * 86400 + (PAXSSNMP - min(PAXSSNMP)) / 80,
      by = .(SEQN, PAXDAYM)]
  # cross-check against the header start time where it is recorded
  hdr <- as.data.table(read_xpt(header_xpt))
  hdr[, SEQN := as.numeric(SEQN)]
  hdr <- hdr[nzchar(trimws(PAXFTIME))]
  hdr[, start_sec := as.numeric(as.POSIXct(paste0("1970-01-01 ", PAXFTIME),
                                           format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))]
  chk <- merge(pax[, .(first_samp = min(PAXSSNMP), d1 = min(as.integer(PAXDAYM))), by = SEQN],
               hdr[, .(SEQN, start_sec)], by = "SEQN")
  # sample 0 was taken at start_sec on day 1, so minute 0 of the first retained day d1 lies
  # (d1 - 1) * 86400 - start_sec seconds after sample 0
  chk[, diff_s := (d1 - 1) * 86400 - start_sec - first_samp / 80]
  cat("  header cross-check (expect ~0 s): max |diff| =", round(max(abs(chk$diff_s)), 1),
      "s over", nrow(chk), "SEQN with a recorded start time\n")

  setorder(pax, SEQN, PAXSSNMP)
  pax[, .(SEQN, PAXDAYM, PAXDAYWM = as.numeric(PAXDAYWM), binary_sleep, UNIX_TIMESTAMP, PAXPREDM)]
}

# ---- run both cycles ------------------------------------------------
cat("== Cycle 2011-12 ==\n")
p12 <- filter_cycle(readRDS(file.path(ddir,"PAXMIN_2011_12.RDS")),
                    file.path(ddir,"PAXHD_G_2011_12.XPT"), "2011-12")
saveRDS(p12, file.path(outdir,"series_1112.rds")); cat("saved series_1112.rds\n\n")
rm(p12); gc()

cat("== Cycle 2013-14 ==\n")
p34 <- filter_cycle(readRDS(file.path(ddir,"PAXMIN_2013_14.RDS")),
                    file.path(ddir,"PAXHD_H_2013_14.XPT"), "2013-14")
saveRDS(p34, file.path(outdir,"series_1314.rds")); cat("saved series_1314.rds\n\n")

# ---- combine + report overlap with analytic sample ------------------
p12 <- readRDS(file.path(outdir,"series_1112.rds"))
ser <- rbind(p12, p34); rm(p12, p34); gc()
saveRDS(ser, file.path(outdir,"series_combined.rds"))
an  <- readRDS(file.path(base,"00_original_paper_materials/data/Combined_data_sets_2_waves_analytic_data_set.RDS"))
cat("filtered SEQN total:", uniqueN(ser$SEQN), "\n")
cat("analytic SEQN total:", length(unique(an$SEQN)), "\n")
cat("overlap with analytic:", length(intersect(ser$SEQN, an$SEQN)), "\n")
cat("DONE part 1\n")
