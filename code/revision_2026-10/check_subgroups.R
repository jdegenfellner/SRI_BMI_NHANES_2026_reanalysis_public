# Day structure of the retained days and the fully adjusted factor (Phillips formula, released series) in two subgroups:
# participants with at least seven consecutive retained days and the participants scored by sleepreg at its default minimum of
# five comparison days. Env: SRI_PROJECT_ROOT, FRAG_DIR (na_raw.rds, na_srE_raw.rds), OUT_DIR (default: working directory).
suppressMessages({library(data.table); library(survey)})
this <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE)
here <- if (is.na(this)) getwd() else dirname(normalizePath(this))
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd()); outdir <- Sys.getenv("OUT_DIR", getwd())
fragdir <- Sys.getenv("FRAG_DIR", file.path(base, "01_recompute_SRI/intermediate/revision_2026-10_results/frag"))
origdata <- file.path(base, "00_original_paper_materials", "data")
d  <- as.data.table(readRDS(file.path(origdata, "df_an1_imputed_MAIN.RDS")))
an <- as.data.table(readRDS(file.path(origdata, "Combined_data_sets_2_waves_analytic_data_set.RDS")))
d  <- merge(d, an[, .(SEQN, WTMEC2YR, SDMVPSU, SDMVSTRA)], by = "SEQN", all.x = TRUE)[, WTMEC4YR := WTMEC2YR / 2]
src <- readLines(file.path(here, "variants4.R")); s0 <- grep("^set.seed\\(3455\\)", src); s1 <- grep("^dm <- merge", src) - 1; eval(parse(text = src[s0:s1]))
ser <- as.data.table(readRDS(file.path(base, "01_recompute_SRI/intermediate/series_combined.rds")))[SEQN %in% d$SEQN, .(SEQN, PAXDAYM)]
days <- unique(ser)[, dayi := as.integer(PAXDAYM)][order(SEQN, dayi)]
nd <- days[, { dd <- dayi; runs <- rle(cumsum(c(1L, diff(dd)) != 1L)); .(ndays = .N, max_consecutive = max(runs$lengths), interior_gap = any(diff(dd) > 1L)) }, by = SEQN]
cat("participants:", nrow(nd), "\n"); cat("share with >= 7 consecutive retained days:", round(mean(nd$max_consecutive >= 7), 3),
    "| with an interior gap:", round(mean(nd$interior_gap), 3), "| fewer than 7 days without interior gap:", round(mean(nd$max_consecutive < 7 & !nd$interior_gap), 3), "\n")
fwrite(nd, file.path(outdir, "participant_ndays.csv"))
g <- merge(readRDS(file.path(fragdir, "na_raw.rds")), readRDS(file.path(fragdir, "na_srE_raw.rds"))[, -"variant"], by = "SEQN", all.x = TRUE); g <- merge(g, nd, by = "SEQN")
res <- list()
for (nm in c("all", "seven_consecutive_days", "sleepreg_default_participants")) {
  sub <- switch(nm, all = g, seven_consecutive_days = g[max_consecutive >= 7], sleepreg_default_participants = g[!is.na(SRI_sleepregE5)])
  dm <- merge(d, sub[, .(SEQN, SRI_hand, SRI_sleepregE5)], by = "SEQN")
  r <- mf_table(dm, dm$SRI_hand, nb = 200, seed = 11); r[, subgroup := nm]; r[, calculator := "Phillips formula"]; res[[nm]] <- r
  if (nm == "sleepreg_default_participants") { r2 <- mf_table(dm, dm$SRI_sleepregE5, nb = 200, seed = 11); r2[, subgroup := nm]; r2[, calculator := "sleepreg, 5 days"]; res[["sleepreg_own"]] <- r2 }
  cat("\n", nm, " n =", nrow(dm), "\n"); print(r[, .(n, M1, M2, M3_ref, M3_other)]) }
out <- rbindlist(res, fill = TRUE); print(out[, .(subgroup, calculator, n, M2, M3_other)]); fwrite(out, file.path(outdir, "check_subgroups.csv"))
