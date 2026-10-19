# =====================================================================
# repro_published_sri.R  (revision 2026-10)
# Reproduce the published SRI values (Degenfellner et al. 2026, JBR) exactly from the
# public PAXMIN minute files, using the two earlier filter states of the original project
# code and the pairing rule of that time (no Saturday->Sunday wrap). Also shows that the
# final archived filter + final pairing rule equals the rebuilt series of this reanalysis.
#
# Inputs:  SRI_PROJECT_ROOT/01_recompute_SRI/data_local/PAXMIN_{2011_12,2013_14}.RDS
#          SRI_PROJECT_ROOT/01_recompute_SRI/intermediate/published_SRI_{2011_12,2013_14}.RDS
#          orig_filters/ (the three states of the filter function of the earlier study, copied verbatim)
# Output:  SRI_PROJECT_ROOT/01_recompute_SRI/intermediate/repro_published_sri.csv
# Runtime: ~10 min, ~12 GB RAM (one cycle in memory at a time).
# =====================================================================
suppressMessages({library(data.table); library(dplyr)})
this <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE)
here <- if (is.na(this)) getwd() else dirname(normalizePath(this))
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd()); filt <- file.path(here, "orig_filters")
ddir <- file.path(base, "01_recompute_SRI", "data_local"); inter <- file.path(base, "01_recompute_SRI", "intermediate")
pubdir <- if (file.exists(file.path(inter, "published_SRI_2011_12.RDS"))) inter else file.path(here, "..", "..", "data")   # the repository ships the published values in data/

# source a filter function verbatim from a file of the original project
source_filter <- function(file) {
  src <- readLines(file); st <- grep("^filter_accelo_data_NHANES_11_to_14 <- function", src)[1]
  en  <- which(src == "}" & seq_along(src) > st)[1]; eval(parse(text = src[st:en]), envir = globalenv()); invisible(c(st, en))
}
# row-shift pairing (minute i vs minute i + 1440 in recording order), as in the original loop
sri_rows <- function(fd, wrap = TRUE) {
  fd <- as.data.table(fd); fd[, wd := as.numeric(PAXDAYWM)]
  fd[, `:=`(nb = shift(binary_sleep, 1440, type = "lead"), nw = shift(wd, 1440, type = "lead")), by = SEQN]
  cons <- if (wrap) fd[, (nw - wd == 1) | (wd == 7 & nw == 1)] else fd[, (nw - wd == 1)]
  fd[, valid := cons & !is.na(binary_sleep) & !is.na(nb)]
  fd[, .(N_v = sum(valid, na.rm = TRUE), delta = sum(valid & binary_sleep == nb, na.rm = TRUE)), by = SEQN
     ][, SRI := fifelse(N_v > 0, -100 + 200 / N_v * delta, NA_real_)][]
}
compare <- function(v, pub, label) {
  m <- merge(v, pub, by = "SEQN", all = TRUE); b <- m[!is.na(SRI) & !is.na(SRI_pub)]; d <- b$SRI - b$SRI_pub
  out <- data.table(variant = label, n = nrow(v), n_published = nrow(pub), common = nrow(b),
                    only_published = sum(is.na(m$SRI) & !is.na(m$SRI_pub)), only_reproduced = sum(!is.na(m$SRI) & is.na(m$SRI_pub)),
                    exact = sum(abs(d) < 1e-6), r = cor(b$SRI, b$SRI_pub))
  print(out); out
}
res <- list()
for (cyc in c("2011_12", "2013_14")) {
  pub <- as.data.table(readRDS(file.path(pubdir, paste0("published_SRI_", cyc, ".RDS"))))[, .(SEQN, SRI_pub = SRI_no_imputation)]
  raw <- readRDS(file.path(ddir, paste0("PAXMIN_", cyc, ".RDS")))
  # (1) final archived filter (R_Codes/2_Read_Clean_Prepare_for_Analysis.R) = Methods as described
  source_filter(file.path(filt, "filter_final.R"))
  fd <- filter_accelo_data_NHANES_11_to_14(raw, "all")$filtered_data
  res[[length(res) + 1]] <- cbind(cycle = cyc, compare(sri_rows(fd, TRUE),  pub, "final filter, final pairing (wrap)"))
  res[[length(res) + 1]] <- cbind(cycle = cyc, compare(sri_rows(fd, FALSE), pub, "final filter, no wrap"))
  rm(fd); gc()
  # (2) earlier filter state that produced the published values of this cycle
  f <- if (cyc == "2011_12") "filter_2011_12_state.R" else "filter_2013_14_state.R"
  source_filter(file.path(filt, f))
  fd <- filter_accelo_data_NHANES_11_to_14(raw, "all")$filtered_data; rm(raw); gc()
  res[[length(res) + 1]] <- cbind(cycle = cyc, compare(sri_rows(fd, TRUE),  pub, paste0(basename(f), ", wrap")))
  res[[length(res) + 1]] <- cbind(cycle = cyc, compare(sri_rows(fd, FALSE), pub, paste0(basename(f), ", no wrap  <- published")))
  rm(fd); gc()
}
out <- rbindlist(res); fwrite(out, file.path(inter, "repro_published_sri.csv")); print(out); cat("DONE\n")
