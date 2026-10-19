# Combines the per-variant model results of na_models_v.R into one table (na_models_all.rds / .csv), the input of fig_forest_grid.R
# and of Table 2 / Supplementary Table S4. Run in the folder that holds na_models_<variant>.rds.
suppressMessages(library(data.table))
vars <- c("raw", "maj1h50", "ggir1h80", "bridge2", "spt2")
tab <- rbindlist(lapply(vars, function(v) readRDS(paste0("na_models_", v, ".rds"))), fill = TRUE)
saveRDS(tab, "na_models_all.rds"); fwrite(tab, "na_models_all.csv"); print(tab[, .(variant, calculator, n, M2)]); cat("DONE\n")
