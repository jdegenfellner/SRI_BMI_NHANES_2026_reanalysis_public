# Model 4 of the earlier study (log BMI on SRI with race/ethnicity x SRI and sex x SRI interactions plus the covariates of model 2,
# survey-weighted) for the published SRI, for every input variant x calculator on the published sample (PART = raw | maj1h50 |
# ggir1h80 | bridge2 | spt2) and for the four series of the random sample (PART = rawsample). For each race/ethnicity x sex group
# the multiplication factor top versus bottom SRI quintile (exp(beta_group * dSRI)) with 95% and 99% percentile intervals from a
# 200-resample participant bootstrap (the published recipe). The full coefficient vectors of every bootstrap fit are stored so that
# predicted BMI across the SRI range (the earlier study's Figure 2) can be derived without refitting.
# Env: SRI_PROJECT_ROOT (working-tree root), PART, NBOOT (default 200), OUT_DIR.
suppressMessages({library(data.table); library(survey)})
base <- Sys.getenv("SRI_PROJECT_ROOT"); part <- Sys.getenv("PART", "published"); nb <- as.integer(Sys.getenv("NBOOT", "200"))
frag <- file.path(base, "01_recompute_SRI/intermediate/revision_2026-10_results/frag"); code <- file.path(base, "SRI_BMI_NHANES_2026_reanalysis/code/revision_2026-10")
rawdir <- file.path(base, "01_recompute_SRI/intermediate/revision_2026-10_results/raw_interim/r1000")
outdir <- Sys.getenv("OUT_DIR", file.path(base, "01_recompute_SRI/intermediate/revision_2026-10_results/model4")); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
origdata <- file.path(base, "00_original_paper_materials", "data")
d  <- as.data.table(readRDS(file.path(origdata, "df_an1_imputed_MAIN.RDS")))
an <- as.data.table(readRDS(file.path(origdata, "Combined_data_sets_2_waves_analytic_data_set.RDS")))
d  <- merge(d, an[, .(SEQN, WTMEC2YR, SDMVPSU, SDMVSTRA)], by = "SEQN", all.x = TRUE)[, WTMEC4YR := WTMEC2YR / 2]
src <- readLines(file.path(code, "variants4.R")); s0 <- grep("^set.seed\\(3455\\)", src); s1 <- grep("^mf_table <- function", src) - 1
prep <- function() eval(parse(text = src[s0:s1]), envir = globalenv())    # set.seed(3455), SRIsim, winz, cov_full: same RNG state as the pipeline
prep()
cov_rest <- setdiff(cov_full, c("RIAGENDR_factor", "RIDRETH3_factor"))
f4 <- as.formula(paste0("log(BMXBMI) ~ SRI*RIDRETH3_factor + SRI*RIAGENDR_factor + ", paste(cov_rest, collapse = " + ")))
eth_lev <- levels(d$RIDRETH3_factor); sex_lev <- levels(d$RIAGENDR_factor)
groups <- CJ(eth = eth_lev, sex = sex_lev, sorted = FALSE)
mkdes <- function(data) svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, nest = TRUE, data = data)
coefs <- function(data) coef(svyglm(f4, design = mkdes(data)))
gbeta <- function(cf) sapply(seq_len(nrow(groups)), function(i) { b <- cf["SRI"]
  for (nm in c(paste0("SRI:RIDRETH3_factor", groups$eth[i]), paste0("RIDRETH3_factor", groups$eth[i], ":SRI"),
               paste0("SRI:RIAGENDR_factor", groups$sex[i]), paste0("RIAGENDR_factor", groups$sex[i], ":SRI"))) if (nm %in% names(cf)) b <- b + cf[nm]
  unname(b) })
m4 <- function(dd, svals, input, calc, seed = 11) {
  dd <- copy(dd); dd[, SRI := winz(svals)]; dd <- dd[is.finite(SRI) & !is.na(BMXBMI) & !is.na(WTMEC4YR)]
  dd[, RIDRETH3_factor := factor(as.character(RIDRETH3_factor), levels = eth_lev)]; dd[, RIAGENDR_factor := factor(as.character(RIAGENDR_factor), levels = sex_lev)]
  med <- dd[, median(SRI), by = .(q = cut(SRI, quantile(SRI, 0:5 / 5), include.lowest = TRUE, labels = 1:5))][order(q), V1]; dSRI <- med[5] - med[1]
  cf0 <- coefs(dd); b0 <- gbeta(cf0); set.seed(seed); n <- nrow(dd)
  CF <- matrix(NA_real_, nb, length(cf0), dimnames = list(NULL, names(cf0)))
  for (b in seq_len(nb)) { bs <- dd[sample(n, replace = TRUE)]; cf <- tryCatch(coefs(bs), error = function(e) NULL)
    if (!is.null(cf)) { k <- intersect(names(cf), colnames(CF)); CF[b, k] <- cf[k] } }
  BB <- apply(CF, 1, gbeta); if (is.null(dim(BB))) BB <- matrix(BB, nrow = nrow(groups))     # groups x nb
  ng <- dd[, .N, by = .(eth = as.character(RIDRETH3_factor), sex = as.character(RIAGENDR_factor))]
  res <- data.table(input = input, calculator = calc, n = n, dSRI = dSRI, medQ1 = med[1], medQ5 = med[5], groups)
  res[, N := ng[res, on = .(eth, sex), N]]
  res[, `:=`(beta = b0, mf = exp(b0 * dSRI),
             lo95 = exp(apply(BB, 1, quantile, .025, na.rm = TRUE) * dSRI), hi95 = exp(apply(BB, 1, quantile, .975, na.rm = TRUE) * dSRI),
             lo99 = exp(apply(BB, 1, quantile, .005, na.rm = TRUE) * dSRI), hi99 = exp(apply(BB, 1, quantile, .995, na.rm = TRUE) * dSRI),
             mf_bootmean = exp(rowMeans(BB, na.rm = TRUE) * dSRI), n_boot_ok = sum(!is.na(CF[, "SRI"])))]
  saveRDS(list(input = input, calculator = calc, n = n, dSRI = dSRI, med = med, coef = cf0, boot_coef = CF, groups = groups),
          file.path(outdir, sprintf("m4_%s_%s.rds", gsub("[^A-Za-z0-9]", "", input), gsub("[^A-Za-z0-9]", "", calc))))
  cat(sprintf("\n[%s / %s] n = %d, dSRI = %.1f, bootstrap fits ok: %d/%d\n", input, calc, n, dSRI, res$n_boot_ok[1], nb))
  print(res[, .(eth, sex, N, mf = round(mf, 3), lo95 = round(lo95, 3), hi95 = round(hi95, 3), lo99 = round(lo99, 3), hi99 = round(hi99, 3))]); flush.console()
  res
}
rows <- list()
if (part == "published") {
  prep(); rows[[1]] <- m4(d, d$SRI_no_imputation_winsorized_rand, "published", "published implementation")
} else if (part == "rawsample") {
  source(file.path(code, "raw_sample.R"))
  res <- readRDS(file.path(rawdir, "raw_derivation_participants.rds")); ord <- fread(file.path(code, "analytic_random_order.tsv"))[, .(SEQN = seqn, rank)]
  excluded <- as.integer(strsplit(Sys.getenv("EXCLUDED_SEQN", "83452,69334,64226"), ",")[[1]])
  samp <- pick_sample(res, ord, "prefix", 1000, excluded)$samp
  sw <- dcast(res[SEQN %in% samp], SEQN ~ deriv, value.var = "sri_A"); dd <- merge(d, sw, by = "SEQN")
  derivs <- c("PAXPREDM", "GGIR_SIB", "GGIR_SPT", "SWaN"); dd <- dd[Reduce(`&`, lapply(derivs, function(s) is.finite(dd[[s]])))]
  cat(sprintf("random sample: %d participants with all four series\n", nrow(dd)))
  for (s in derivs) { prep(); rows[[s]] <- m4(dd, dd[[s]], "random sample", s) }
} else {
  g <- merge(readRDS(file.path(frag, paste0("na_", part, ".rds"))), readRDS(file.path(frag, paste0("na_srE_", part, ".rds")))[, -"variant"], by = "SEQN", all.x = TRUE)
  dm <- merge(d, g, by = "SEQN")
  for (cc in c("SRI_hand", "SRI_GGIR", "SRI_sleepregE5", "SRI_sleepregE3")) { prep(); rows[[cc]] <- m4(dm, dm[[cc]], part, cc) }
}
fwrite(rbindlist(rows), file.path(outdir, paste0("m4_", part, ".csv"))); cat("\nDONE", part, "\n")
