# swan_prep.R <extract_dir> <out_dir> : NHANES-Stunden-CSV -> mHealth-Layout (<out>/<YYYY-MM-DD>/<HH>/...sensor.csv), Werte in g, ohne ISM-Umkehr (macht keinen Unterschied, Pilot 62164)
suppressMessages(library(data.table)); setDTthreads(2)
a <- commandArgs(trailingOnly = TRUE); src <- a[1]; out <- a[2]
fl <- sort(list.files(src, pattern = "sensor\\.csv$", full.names = TRUE, recursive = TRUE)); fl <- fl[!grepl("/\\._", fl)]
for (f in fl) {
  x <- fread(f, colClasses = c("character", "numeric", "numeric", "numeric"))
  setnames(x, c("HEADER_TIME_STAMP", "X_ACCELERATION_METERS_PER_SECOND_SQUARED", "Y_ACCELERATION_METERS_PER_SECOND_SQUARED", "Z_ACCELERATION_METERS_PER_SECOND_SQUARED"))
  d <- substr(x$HEADER_TIME_STAMP[1], 1, 10); h <- substr(x$HEADER_TIME_STAMP[1], 12, 13)
  od <- file.path(out, d, h); dir.create(od, recursive = TRUE, showWarnings = FALSE)
  x[, HEADER_TIME_STAMP := ifelse(nchar(HEADER_TIME_STAMP) == 19, paste0(HEADER_TIME_STAMP, ".000"), HEADER_TIME_STAMP)]
  fwrite(x, file.path(od, sprintf("AndroidWearWatch-AccelerationCalibrated-NA.%s-%s-00-00-000-P0000.sensor.csv", d, h)))
}
cat("Stunden:", length(fl), "\n")
