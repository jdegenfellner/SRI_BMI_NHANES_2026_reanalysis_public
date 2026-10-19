# Slim copy of PAXMIN for the raw-signal comparison: analytic participants only, minutes with PAXTSM == 60 & PAXQFM == 0 (as in raw_derivation_participants.R), needed columns only.
suppressMessages(library(data.table))
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())   # working-tree root (contains 00_original_paper_materials/ and 01_recompute_SRI/)
an <- unique(as.data.table(readRDS(file.path(base, "00_original_paper_materials/data/df_an1_imputed_MAIN.RDS")))$SEQN)
for (cyc in c("2011_12", "2013_14")) {
  t0 <- Sys.time()
  x <- readRDS(file.path(base, sprintf("01_recompute_SRI/data_local/PAXMIN_%s.RDS", cyc))); setDT(x)
  cat(cyc, "loaded:", nrow(x), "Zeilen,", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
  x <- x[SEQN %in% an & PAXTSM == 60 & PAXQFM == 0, .(SEQN = as.integer(SEQN), PAXDAYM = as.integer(PAXDAYM), PAXSSNMP, PAXTSM, PAXQFM, PAXAISMM, PAXPREDM = as.character(PAXPREDM))]
  gc(verbose = FALSE)
  saveRDS(x, file.path(base, sprintf("01_recompute_SRI/intermediate/slim_paxmin/PAXMIN_slim_%s.rds", cyc)))
  cat(cyc, "slim:", nrow(x), "rows,", uniqueN(x$SEQN), "SEQN,", round(object.size(x) / 2^20), "MB,", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
  rm(x); gc(verbose = FALSE)
}
