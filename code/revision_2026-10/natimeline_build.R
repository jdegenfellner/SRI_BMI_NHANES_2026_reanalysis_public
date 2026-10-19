# Build the faithful NA timeline: one row per minute from monitor placement to the end of the last retained day,
# with the retained days' states and NA (invalid) elsewhere. Timestamps = seconds since 00:00 of the placement day.
suppressMessages(library(data.table))
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())   # working tree root (contains 00_original_paper_materials/ and 01_recompute_SRI/)
an <- unique(as.data.table(readRDS(file.path(base, "00_original_paper_materials/data/df_an1_imputed_MAIN.RDS")))$SEQN)
ser <- as.data.table(readRDS(file.path(base, "01_recompute_SRI/intermediate/series_combined.rds")))[SEQN %in% an]
setorder(ser, SEQN, UNIX_TIMESTAMP)
ser[, dayi := as.integer(PAXDAYM)]; ser[, day_start := UNIX_TIMESTAMP[1], by = .(SEQN, dayi)]
ser[, ts := (dayi - 1) * 86400 + round(UNIX_TIMESTAMP - day_start)]
# placement time: old UNIX_TIMESTAMP was PAXSSNMP/80 (seconds since placement), so placement = ts - UNIX_TIMESTAMP at the first retained minute
pl <- ser[, .(placement = ts[1] - UNIX_TIMESTAMP[1], ts_end = max(ts) + 60, d1 = min(dayi), nd = uniqueN(dayi)), by = SEQN]
pl[, placement := 60 * floor(placement / 60)]
cat("placement time-of-day (h):\n"); print(round(quantile((pl$placement %% 86400)/3600, c(0, .05, .25, .5, .75, .95, 1)), 2))
cat("participants placed exactly at midnight:", sum(pl$placement %% 86400 == 0), "\n")
pl[placement %% 86400 == 0, placement := placement + 30]        # sleepreg cannot start a recording at exactly 00:00:00; shift those few by 30 s
# full minute grid per participant
grid <- pl[, .(ts = seq(placement, ts_end - 60, by = 60)), by = SEQN]
full <- merge(grid, ser[, .(SEQN, ts, binary_sleep, PAXDAYWM, dayi)], by = c("SEQN", "ts"), all.x = TRUE)
setorder(full, SEQN, ts)
cat("rows:", nrow(full), " retained minutes:", sum(!is.na(full$dayi)), " participants:", uniqueN(full$SEQN), "\n")
saveRDS(full, "series_natimeline.rds"); saveRDS(pl, "placement.rds")
# attach the smoothed input variants (keyed by SEQN, ts_fix == ts)
v <- as.data.table(readRDS("series_variants.rds"))[, .(SEQN, ts = ts_fix, maj1h50, ggir1h80, bridge2, spt2)]
full <- merge(full, v, by = c("SEQN", "ts"), all.x = TRUE); setorder(full, SEQN, ts)
stopifnot(all(is.na(full$maj1h50) | !is.na(full$dayi)))
saveRDS(full, "series_natimeline_variants.rds"); cat("saved series_natimeline_variants.rds\n")
