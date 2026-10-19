# ggir_one.R <seqn> <pax_g|pax_h> : NHANES-Rohdaten eines Teilnehmers -> GGIR 1-5 (HDCZA, vanHees2015 SIB, ZC-Metriken fuer spaetere Varianten)
# behaelt nur: part2/part4/part5-Zusammenfassungen, Epochenserie (ms5.outraw), QC; loescht die entpackte CSV.
suppressMessages({library(GGIR); library(data.table)}); data.table::setDTthreads(2)
args <- commandArgs(trailingOnly = TRUE); seqn <- args[1]; cyc <- args[2]
base <- Sys.getenv("NH_BASE", "/Volumes/LaCie/NHANES_1112_1314"); out <- file.path(base, "_ggir", cyc); dir.create(out, showWarnings = FALSE, recursive = TRUE)
done <- file.path(out, "done", paste0(seqn, ".txt")); dir.create(dirname(done), showWarnings = FALSE)
if (file.exists(done)) { cat("schon fertig\n"); quit(save = "no") }
raw <- Sys.getenv("NH_RAW_FILE", file.path(base, if (cyc == "pax_g") "raw_accel_data_201112" else "raw_accel_data_201314", paste0(seqn, ".tar.bz2")))
W <- file.path(Sys.getenv("NH_TMP", file.path(Sys.getenv("HOME"), ".ggir_tmp")), paste0("ggir_", seqn)); ex <- file.path(W, "extract"); acc <- file.path(W, "acc")   # interne SSD statt ExFAT
unlink(W, recursive = TRUE); dir.create(ex, recursive = TRUE); dir.create(acc)
t0 <- Sys.time()
system2("tar", c("-xjf", shQuote(raw), "-C", shQuote(ex)))
fl <- sort(list.files(ex, pattern = "sensor\\.csv$", full.names = TRUE, recursive = TRUE))
dt <- rbindlist(lapply(fl, fread, colClasses = c("character", "numeric", "numeric", "numeric")))
fwrite(dt, file.path(acc, paste0(seqn, ".csv"))); start_ts <- dt$HEADER_TIMESTAMP[1]; n_samp <- nrow(dt); rm(dt); unlink(ex, recursive = TRUE)
t1 <- Sys.time()
if (nzchar(Sys.getenv("NH_GGIR_PATCH"))) source(Sys.getenv("NH_GGIR_PATCH"))
GGIR(mode = 1:5, datadir = acc, outputdir = W, studyname = paste0("nhanes_", seqn),
     rmc.firstrow.acc = 2, rmc.col.acc = 2:4, rmc.col.time = 1, rmc.unit.time = "POSIX",
     rmc.format.time = "%Y-%m-%d %H:%M:%OS", rmc.unit.acc = "g", rmc.sf = 80, desiredtz = "UTC", configtz = "UTC",
     windowsizes = c(5, 900, 3600), do.cal = TRUE, do.enmo = TRUE, do.anglez = TRUE, do.zcx = TRUE, do.zcy = TRUE, do.zcz = TRUE,
     HASPT.algo = "HDCZA", HASIB.algo = "vanHees2015", def.noc.sleep = 1, sleepwindowType = "SPT",
     outliers.only = FALSE, includenightcrit = 16, includedaycrit = 16,
     do.sibreport = TRUE, save_ms5rawlevels = TRUE, save_ms5raw_format = "csv", save_ms5raw_without_invalid = FALSE,
     do.report = c(2, 4, 5), visualreport = FALSE, overwrite = TRUE, do.parallel = FALSE, verbose = FALSE)
t2 <- Sys.time()
# Ergebnisse einsammeln
og <- file.path(W, "output_acc"); keep <- file.path(out, seqn); dir.create(keep, showWarnings = FALSE)
for (f in c("results/part2_daysummary.csv", "results/part4_nightsummary_sleep_cleaned.csv", "results/QC/part4_nightsummary_sleep_full.csv",
            "results/part5_daysummary_MM_L40M100V400_T5A5.csv", "results/QC/data_quality_report.csv")) if (file.exists(file.path(og, f))) file.copy(file.path(og, f), file.path(keep, basename(f)), overwrite = TRUE)
ts <- list.files(file.path(og, "meta", "ms5.outraw"), pattern = "\\.csv$", full.names = TRUE, recursive = TRUE); ts <- ts[!grepl("behavioralcodes", ts)]
if (length(ts)) { x <- fread(ts[1]); fwrite(x, file.path(keep, paste0(seqn, "_epochs.csv.gz"))) }
bc <- list.files(file.path(og, "meta", "ms5.outraw"), pattern = "behavioralcodes", full.names = TRUE); if (length(bc)) file.copy(bc[1], file.path(keep, "behavioralcodes.csv"), overwrite = TRUE)
file.copy(list.files(file.path(og, "meta", "ms3.out"), full.names = TRUE)[1], file.path(keep, paste0(seqn, "_ms3.RData")), overwrite = TRUE)  # SIB-Daten fuer spaetere Varianten
file.copy(list.files(file.path(og, "meta", "basic"), full.names = TRUE)[1], file.path(keep, paste0(seqn, "_basic.RData")), overwrite = TRUE) # Part-1-Metriken (inkl. ZC) fuer Varianten
writeLines(c(paste0("seqn\t", seqn), paste0("cycle\t", cyc), paste0("raw_start\t", start_ts), paste0("n_samples\t", n_samp),
             paste0("prep_min\t", round(as.numeric(difftime(t1, t0, units = "mins")), 2)), paste0("ggir_min\t", round(as.numeric(difftime(t2, t1, units = "mins")), 2))), done)
unlink(W, recursive = TRUE)
cat(seqn, "fertig: prep", round(as.numeric(difftime(t1, t0, units = "mins")), 1), "min, GGIR", round(as.numeric(difftime(t2, t1, units = "mins")), 1), "min\n")
