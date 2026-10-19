# Hand / GGIR (frac_valid-filtered) / sleepreg on the NA timeline for ONE input variant. Usage: Rscript na_grid.R <raw|maj1h50|ggir1h80|bridge2|spt2>
suppressMessages({library(data.table); library(GGIR)})
v <- commandArgs(trailingOnly = TRUE)[1]; stopifnot(v %in% c("raw","maj1h50","ggir1h80","bridge2","spt2"))
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())   # working tree root (contains 00_original_paper_materials/ and 01_recompute_SRI/)
source(file.path(base, "01_recompute_SRI/sleepreg_src/sleepreg_allFunctions.R")); load(file.path(base, "01_recompute_SRI/sleepreg_src/sysdata.rda"))
here <- getwd(); col <- if (v == "raw") "binary_sleep" else v
ser <- as.data.table(readRDS("series_natimeline_variants.rds"))[, c("SEQN", "ts", col), with = FALSE]; setnames(ser, col, "s")
setorder(ser, SEQN, ts); setkey(ser, SEQN); seqns <- sort(unique(ser$SEQN))
hand_na <- function(s, M = 1440L) { n <- length(s); i <- 1:(n - M); a <- s[i]; b <- s[i + M]; ok <- !is.na(a) & !is.na(b); if (!any(ok)) return(NA_real_); -100 + 200 * mean(a[ok] == b[ok]) }
ggir_na <- function(ts, s) { df <- data.frame(time = as.POSIXct(ts, origin = "1970-01-01", tz = "UTC"), sleepstate = as.numeric(s), invalid = as.integer(is.na(s)))
  out <- tryCatch(GGIR::CalcSleepRegularityIndex(data = df, epochsize = 60, desiredtz = "UTC"), error = function(e) NA)
  if (!is.data.frame(out)) return(c(filt = NA, all = NA, npairs = NA))
  ok <- out$frac_valid >= 0.66
  c(filt = if (any(ok)) mean(out$SleepRegularityIndex[ok]) else NA_real_, all = mean(out$SleepRegularityIndex, na.rm = TRUE), npairs = sum(ok)) }
res <- data.table(SEQN = seqns, SRI_hand = NA_real_, SRI_GGIR = NA_real_, SRI_GGIR_unfiltered = NA_real_, GGIR_pairs = NA_real_); t0 <- Sys.time()
for (k in seq_along(seqns)) { x <- ser[.(seqns[k])]; res$SRI_hand[k] <- hand_na(x$s); g <- ggir_na(x$ts, x$s); res$SRI_GGIR[k] <- g["filt"]; res$SRI_GGIR_unfiltered[k] <- g["all"]; res$GGIR_pairs[k] <- g["npairs"]
  if (k %% 1000 == 0) cat(sprintf("  [%s] hand/GGIR %d/%d (%.1f min)\n", v, k, length(seqns), as.numeric(difftime(Sys.time(), t0, units = "mins")))) }
saveRDS(res, paste0("na_", v, "_hand_ggir.rds")); cat("[", v, "] hand+GGIR done\n")
bindir <- file.path(here, paste0("na_sr_in_", v)); unlink(bindir, recursive = TRUE); dir.create(bindir)
for (q in seqns) { x <- ser[.(q)]; fwrite(data.table(binary_sleep = x$s, ts = x$ts), file.path(bindir, paste0("sw_", q, ".csv"))) }
invisible(capture.output(SWS_from_binarySW(binarySWdir = bindir, tsCol = 2, binaryCol = 1)))
for (md in c(5, 3)) { od <- file.path(here, paste0("na_sr_out_", v, "_min", md)); unlink(od, recursive = TRUE); dir.create(od)
  invisible(capture.output(SRI_from_binary(binarydir = file.path(bindir, "SWS_output"), alloutdir = od, col.trans = 1, col.timestamp = 2, wr.raster = FALSE, overwr = TRUE, minSRIdays = md)))
  sr <- fread(file.path(od, "SRI.csv")); sr[, SEQN := as.numeric(gsub("\\D", "", File))]; sr <- sr[, .(SEQN, SRI = as.numeric(SRI), days = as.numeric(SRI_days))]
  setnames(sr, c("SRI", "days"), paste0(c("SRI_sleepreg", "sleepreg_days"), md)); res <- merge(res, sr, by = "SEQN", all.x = TRUE) }
res[, variant := v]; saveRDS(res, paste0("na_", v, ".rds")); unlink(bindir, recursive = TRUE)
cat(sprintf("[%s] DONE: hand median %.1f | GGIR(filt) %.1f | GGIR(unfilt) %.1f | sleepreg5 %.1f (n=%d) | sleepreg3 %.1f (n=%d) | r(hand,GGIR) %.3f r(hand,sleepreg5) %.3f\n", v,
  median(res$SRI_hand, na.rm=T), median(res$SRI_GGIR, na.rm=T), median(res$SRI_GGIR_unfiltered, na.rm=T), median(res$SRI_sleepreg5, na.rm=T), sum(!is.na(res$SRI_sleepreg5)), median(res$SRI_sleepreg3, na.rm=T), sum(!is.na(res$SRI_sleepreg3)),
  cor(res$SRI_hand, res$SRI_GGIR, use="pair"), cor(res$SRI_hand, res$SRI_sleepreg5, use="pair")))
