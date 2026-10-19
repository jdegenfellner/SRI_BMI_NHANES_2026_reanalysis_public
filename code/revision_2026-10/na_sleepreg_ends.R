# sleepreg with the recording ENDS ABSENT (timeline from 00:00 of the first retained day to 24:00 of the last, interior gaps NA),
# timestamps shifted +30 s for sleepreg's origin-clock parser. Usage: Rscript na_sleepreg_ends.R <variant>
suppressMessages(library(data.table))
v <- commandArgs(trailingOnly = TRUE)[1]; stopifnot(v %in% c("raw","maj1h50","ggir1h80","bridge2","spt2"))
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())   # working tree root (contains 00_original_paper_materials/ and 01_recompute_SRI/)
source(file.path(base, "01_recompute_SRI/sleepreg_src/sleepreg_allFunctions.R")); load(file.path(base, "01_recompute_SRI/sleepreg_src/sysdata.rda"))
here <- getwd(); col <- if (v == "raw") "binary_sleep" else v
ser <- as.data.table(readRDS("series_natimeline_variants.rds"))[, c("SEQN", "ts", "dayi", col), with = FALSE]; setnames(ser, col, "s")
# trim to [00:00 of first retained day, 24:00 of last retained day]
lim <- ser[!is.na(dayi), .(t0 = (min(dayi) - 1) * 86400, t1 = max(dayi) * 86400), by = SEQN]
ser <- merge(ser, lim, by = "SEQN")[ts >= t0 & ts < t1][, c("t0","t1","dayi") := NULL]; setorder(ser, SEQN, ts); setkey(ser, SEQN); seqns <- sort(unique(ser$SEQN))
bindir <- file.path(here, paste0("na_srE_in_", v)); unlink(bindir, recursive = TRUE); dir.create(bindir)
for (q in seqns) { x <- ser[.(q)]; fwrite(data.table(binary_sleep = x$s, ts = x$ts + 30), file.path(bindir, paste0("sw_", q, ".csv"))) }
invisible(capture.output(SWS_from_binarySW(binarySWdir = bindir, tsCol = 2, binaryCol = 1)))
res <- data.table(SEQN = seqns)
for (md in c(5, 3)) { od <- file.path(here, paste0("na_srE_out_", v, "_min", md)); unlink(od, recursive = TRUE); dir.create(od)
  invisible(capture.output(SRI_from_binary(binarydir = file.path(bindir, "SWS_output"), alloutdir = od, col.trans = 1, col.timestamp = 2, wr.raster = FALSE, overwr = TRUE, minSRIdays = md)))
  sr <- fread(file.path(od, "SRI.csv")); sr[, SEQN := as.numeric(gsub("\\D", "", File))]; sr <- sr[, .(SEQN, SRI = as.numeric(SRI), days = as.numeric(SRI_days))]
  setnames(sr, c("SRI", "days"), paste0(c("SRI_sleepregE", "sleepregE_days"), md)); res <- merge(res, sr, by = "SEQN", all.x = TRUE) }
res[, variant := v]; saveRDS(res, paste0("na_srE_", v, ".rds")); unlink(bindir, recursive = TRUE)
cat(sprintf("[%s] DONE ends-absent: sleepreg5 median %.1f (n=%d) | sleepreg3 median %.1f (n=%d)\n", v, median(res$SRI_sleepregE5, na.rm=T), sum(!is.na(res$SRI_sleepregE5)), median(res$SRI_sleepregE3, na.rm=T), sum(!is.na(res$SRI_sleepregE3))))
