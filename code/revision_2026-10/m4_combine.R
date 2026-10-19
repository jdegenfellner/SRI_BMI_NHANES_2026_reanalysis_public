# Combines the per-input Model 4 results of na_models_m4.R (m4_<part>.csv) into m4_all.csv, the input of m4_tables.R, and prints
# the group factors against the earlier study's Table 3. Env: SRI_PROJECT_ROOT (working-tree root); M4_DIR overrides the results folder.
suppressMessages(library(data.table))
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())
dir <- Sys.getenv("M4_DIR", file.path(base, "01_recompute_SRI/intermediate/revision_2026-10_results/model4"))
files <- setdiff(list.files(dir, pattern = "^m4_.*[.]csv$", full.names = TRUE), file.path(dir, c("m4_all.csv", "m4_predictions.csv")))
x <- rbindlist(lapply(files, fread), fill = TRUE); x <- x[input != ""]; fwrite(x, file.path(dir, "m4_all.csv"))
x[, cell := sprintf("%.3f (%.3f, %.3f)", mf, lo95, hi95)]
ethshort <- c("Mexican American" = "MexAm", "Other Hispanic" = "OthHisp", "Non-Hispanic White" = "NHWhite", "Non-Hispanic Black" = "NHBlack", "Non-Hispanic Asian" = "NHAsian", "Other Race - Including Multi-Racial" = "Other")
x[, g := paste(ethshort[eth], substr(sex, 1, 1))]
cat("== Published SRI, refit (earlier study, Table 3, Model 4, Q1->Q5: MexAm M 1.04 F 0.98 | OthHisp M 0.991 F 0.934 | NHWhite M 0.967 F 0.91 | NHBlack M 1.036 F 0.975 | NHAsian M 0.979 F 0.922 | Other M 0.937 F 0.882)\n")
print(x[input == "published", .(g, N, cell)], nrows = 20)
cat("\n== Phillips formula by input (published sample)\n")
print(dcast(x[calculator %in% c("published implementation", "SRI_hand")], g ~ factor(input, levels = c("published", "raw", "maj1h50", "ggir1h80", "bridge2", "spt2")), value.var = "cell"), nrows = 20)
cat("\n== Random sample: factor by series\n")
print(dcast(x[input == "random sample"], g + N ~ factor(calculator, levels = c("PAXPREDM", "GGIR_SIB", "GGIR_SPT", "SWaN")), value.var = "cell"), nrows = 20)
cat("written:", file.path(dir, "m4_all.csv"), "\n")
