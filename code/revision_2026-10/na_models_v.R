# Table-3 recipe models for ONE input variant x {hand, GGIR, sleepregE5, sleepregE3}; 200 bootstrap. Usage: Rscript na_models_v.R <variant>
suppressMessages({library(data.table); library(survey)})
v <- commandArgs(trailingOnly = TRUE)[1]
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())   # working tree root (contains 00_original_paper_materials/ and 01_recompute_SRI/)
origdata <- file.path(base, "00_original_paper_materials", "data")
d  <- as.data.table(readRDS(file.path(origdata, "df_an1_imputed_MAIN.RDS")))
an <- as.data.table(readRDS(file.path(origdata, "Combined_data_sets_2_waves_analytic_data_set.RDS")))
d  <- merge(d, an[, .(SEQN, WTMEC2YR, SDMVPSU, SDMVSTRA)], by = "SEQN", all.x = TRUE)[, WTMEC4YR := WTMEC2YR/2]
g <- merge(readRDS(paste0("na_", v, ".rds")), readRDS(paste0("na_srE_", v, ".rds"))[, -"variant"], by = "SEQN", all.x = TRUE)
this <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE); here <- if (is.na(this)) getwd() else dirname(normalizePath(this)); src <- readLines(file.path(here, "variants4.R")); s0 <- grep("^set.seed\\(3455\\)", src); s1 <- grep("^dm <- merge", src) - 1; eval(parse(text = src[s0:s1]))
dm <- merge(d, g, by = "SEQN"); rows <- list()
for (cc in c("SRI_hand","SRI_GGIR","SRI_sleepregE5","SRI_sleepregE3","SRI_sleepreg5","SRI_sleepreg3")) { rows[[cc]] <- cbind(variant = v, calculator = cc, mf_table(dm, dm[[cc]])); cat("[", v, "]", cc, "done\n") }
tab <- rbindlist(rows); saveRDS(tab, paste0("na_models_", v, ".rds")); print(tab[, .(calculator, n, dSRI, M2, M3_other)]); cat("DONE\n")
