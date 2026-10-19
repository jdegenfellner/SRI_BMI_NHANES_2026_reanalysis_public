# Four implementations on identical NA-free input: the 300 NHANES participants with seven consecutive retained days, a mean sleep share
# between 0.15 and 0.6 and the fewest missing minutes (non-wear or unknown; at most 2.8% of minutes), filled by last observation carried
# forward and, at the start, backward (identical for all four). sleepreg: ends absent, +30 s.
suppressMessages({library(data.table); library(GGIR); library(zoo)})
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())   # working tree root (contains 00_original_paper_materials/ and 01_recompute_SRI/)
source(file.path(base, "01_recompute_SRI/sleepreg_src/sleepreg_allFunctions.R")); load(file.path(base, "01_recompute_SRI/sleepreg_src/sysdata.rda"))
s <- as.data.table(readRDS("series_clockfixed.rds"))
ok <- s[, .(n = .N, nna = sum(is.na(binary_sleep)), nd = uniqueN(dayi), consec = all(diff(sort(unique(dayi))) == 1), ms = mean(binary_sleep, na.rm=TRUE)), by = SEQN][n == 7*1440 & nd == 7 & consec & ms > 0.15 & ms < 0.6][order(nna)]
cat("7-day participants:", nrow(ok), "; NA minutes quantiles:", quantile(ok$nna, c(0,.1,.5,.9)), "\n"); pick <- head(ok$SEQN, 300); cat("300 lowest-NA participants: max NA minutes", max(head(ok$nna, 300)), "= ", round(100*max(head(ok$nna,300))/10080, 2), "% of minutes\n")
x <- s[SEQN %in% pick, .(SEQN, ts = ts_fix, s = binary_sleep)][order(SEQN, ts)]; rm(s); gc()
x[, s := na.locf(na.locf(s, na.rm = FALSE), fromLast = TRUE), by = SEQN]; stopifnot(!anyNA(x$s))

fwrite(x, "xcheck4_series.csv")
hand <- function(v, M = 1440L) { n <- length(v); i <- 1:(n-M); -100 + 200*mean(v[i] == v[i+M]) }
ggir <- function(ts, v) { df <- data.frame(time = as.POSIXct(ts, origin="1970-01-01", tz="UTC"), sleepstate = as.numeric(v), invalid = 0L); o <- GGIR::CalcSleepRegularityIndex(data = df, epochsize = 60, desiredtz = "UTC"); mean(o$SleepRegularityIndex[o$frac_valid >= 0.66]) }
r <- x[, .(SRI_hand = hand(s), SRI_GGIR = ggir(ts, s)), by = SEQN]
bd <- file.path(getwd(), "xcheck4_sr"); unlink(bd, recursive = TRUE); dir.create(bd)
for (q in pick) { y <- x[SEQN == q]; fwrite(data.table(b = y$s, t = y$ts + 30), file.path(bd, paste0("sw_", q, ".csv"))) }
invisible(capture.output({ SWS_from_binarySW(binarySWdir = bd, tsCol = 2, binaryCol = 1); od <- file.path(bd, "out"); dir.create(od); SRI_from_binary(binarydir = file.path(bd, "SWS_output"), alloutdir = od, col.trans = 1, col.timestamp = 2, wr.raster = FALSE, overwr = TRUE, minSRIdays = 5) }))
sr <- fread(file.path(od, "SRI.csv")); sr[, SEQN := as.numeric(gsub("\\D", "", File))]; r <- merge(r, sr[, .(SEQN, SRI_sleepreg = as.numeric(SRI))], by = "SEQN")
saveRDS(r, "xcheck4_r.rds"); cat("R side done, n =", nrow(r), "\n")
