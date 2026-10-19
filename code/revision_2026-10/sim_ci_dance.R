# "Dance of the confidence intervals": take the fitted fully adjusted model (Model 2) or the sex-interaction model (Model 3)
# on the released series as the true data-generating process, regenerate log-BMI for the same participants (fitted value plus a
# residual drawn from the empirical residuals), refit the model with the published recipe (survey-weighted svyglm, multiplication
# factor top versus bottom SRI quintile, 200-resample bootstrap) and record the 95% and 99% percentile intervals. Repeated NSIM times.
# Env: SRI_PROJECT_ROOT (working-tree root), MODEL = M2 | M3w (women), NSIM (default 20), NBOOT (default 200), OUT_DIR.
suppressMessages({library(data.table); library(survey)})
base <- Sys.getenv("SRI_PROJECT_ROOT"); model <- Sys.getenv("MODEL", "M2"); nsim <- as.integer(Sys.getenv("NSIM", "20")); nb <- as.integer(Sys.getenv("NBOOT", "200"))
frag <- file.path(base, "01_recompute_SRI/intermediate/revision_2026-10_results/frag"); code <- file.path(base, "SRI_BMI_NHANES_2026_reanalysis/code/revision_2026-10")
outdir <- Sys.getenv("OUT_DIR", file.path(base, "01_recompute_SRI/intermediate/revision_2026-10_results/sim_ci_dance")); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
origdata <- file.path(base, "00_original_paper_materials", "data")
d  <- as.data.table(readRDS(file.path(origdata, "df_an1_imputed_MAIN.RDS")))
an <- as.data.table(readRDS(file.path(origdata, "Combined_data_sets_2_waves_analytic_data_set.RDS")))
d  <- merge(d, an[, .(SEQN, WTMEC2YR, SDMVPSU, SDMVSTRA)], by = "SEQN", all.x = TRUE)[, WTMEC4YR := WTMEC2YR / 2]
src <- readLines(file.path(code, "variants4.R")); s0 <- grep("^set.seed\\(3455\\)", src); s1 <- grep("^mf_table <- function", src) - 1
eval(parse(text = src[s0:s1]))                      # set.seed(3455), SRIsim, winz, cov_full, cov_full_nosex, fit_set (same RNG state as the pipeline)
dm <- merge(d, readRDS(file.path(frag, "na_raw.rds")), by = "SEQN")
dd <- copy(dm); dd[, SRI := winz(SRI_hand)]; dd <- dd[is.finite(SRI) & !is.na(BMXBMI) & !is.na(WTMEC4YR)]
med <- dd[, median(SRI), by = .(q = cut(SRI, quantile(SRI, 0:5 / 5), include.lowest = TRUE, labels = 1:5))][order(q), V1]; dSRI <- med[5] - med[1]
f2 <- as.formula(paste0("log(BMXBMI) ~ SRI + ", paste(cov_full, collapse = " + ")))
f3 <- as.formula(paste0("log(BMXBMI) ~ RIAGENDR_factor*SRI - SRI + SRI + ", paste(cov_full_nosex, collapse = " + ")))
mkdes <- function(data) svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, nest = TRUE, data = data)
beta <- function(data) { cf <- coef(svyglm(if (model == "M2") f2 else f3, design = mkdes(data)))
  if (model == "M2") return(unname(cf["SRI"]))
  ix <- grep("RIAGENDR_factor.*:SRI|SRI:RIAGENDR_factor", names(cf), value = TRUE); unname(cf["SRI"] + cf[ix[1]]) }   # women = second level (Female)
truth_fit <- svyglm(if (model == "M2") f2 else f3, design = mkdes(dd))
truth <- exp(beta(dd) * dSRI); mu <- fitted(truth_fit); e <- residuals(truth_fit, type = "response")
cat(sprintf("model %s: n = %d, dSRI = %.1f, true MF = %.4f, residual SD = %.3f, levels %s\n", model, nrow(dd), dSRI, truth, sd(e), paste(levels(dd$RIAGENDR_factor), collapse = "/")))
out <- file.path(outdir, paste0("sim_", model, ".csv")); if (file.exists(out)) file.remove(out)
for (s in seq_len(nsim)) {
  set.seed(1000 + s); ds <- copy(dd); ds[, BMXBMI := exp(mu + sample(e, replace = TRUE))]
  pt <- exp(beta(ds) * dSRI); n <- nrow(ds)
  B <- replicate(nb, { bs <- ds[sample(n, replace = TRUE)]; tryCatch(exp(beta(bs) * dSRI), error = function(err) NA_real_) })
  r <- data.table(model = model, sim = s, truth = truth, mf = pt, lo95 = quantile(B, .025, na.rm = TRUE), hi95 = quantile(B, .975, na.rm = TRUE),
                  lo99 = quantile(B, .005, na.rm = TRUE), hi99 = quantile(B, .995, na.rm = TRUE), sd_log = sd(log(B), na.rm = TRUE), n_ok = sum(!is.na(B)))
  fwrite(r, out, append = file.exists(out)); cat(sprintf("sim %2d/%d: MF %.3f  95%% (%.3f, %.3f)  99%% (%.3f, %.3f)\n", s, nsim, pt, r$lo95, r$hi95, r$lo99, r$hi99)); flush.console()
}
cat("DONE\n")
