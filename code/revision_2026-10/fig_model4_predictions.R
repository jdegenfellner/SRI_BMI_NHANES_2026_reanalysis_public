# Predicted BMI across the SRI range from model 4 (the earlier study's Figure 2), for every input/calculator stored by na_models_m4.R.
# Covariates are held at the published sample's means (numeric) and most frequent levels (factors), as in the earlier study; the
# 95% band is the percentile range of the predictions over the stored bootstrap coefficient vectors. Writes a prediction table and
# a multipage PDF with one page per input/calculator plus PNGs of selected pages.
# Env: SRI_PROJECT_ROOT; optional ONLY (comma-separated list of m4_*.rds basenames to plot).
suppressMessages(library(data.table))
base <- Sys.getenv("SRI_PROJECT_ROOT"); m4dir <- file.path(base, "01_recompute_SRI/intermediate/revision_2026-10_results/model4")
origdata <- file.path(base, "00_original_paper_materials", "data")
d <- as.data.table(readRDS(file.path(origdata, "df_an1_imputed_MAIN.RDS")))
cov_rest <- c("RIDAGEYR", "DMDEDUC2_factor", "INDHHIN2_factor", "OCD150_factor", "alcohol_consumption", "smoking_status", "LBXVIDMS", "DR1TKCAL_winsorized", "Depression_score_PHQ_9", "DMDMARTL_factor", "activity_level")
fixed <- lapply(cov_rest, function(v) { x <- d[[v]]; if (is.factor(x)) factor(names(sort(table(x), decreasing = TRUE))[1], levels = levels(x)) else mean(x, na.rm = TRUE) }); names(fixed) <- cov_rest
eth_lev <- levels(d$RIDRETH3_factor); sex_lev <- levels(d$RIAGENDR_factor)
grid <- CJ(SRI = seq(0, 100, by = 2.5), eth = eth_lev, sex = sex_lev, sorted = FALSE)
nd <- cbind(grid, as.data.table(fixed)[rep(1, nrow(grid))]); nd[, RIDRETH3_factor := factor(eth, levels = eth_lev)]; nd[, RIAGENDR_factor := factor(sex, levels = sex_lev)]
X <- model.matrix(as.formula(paste0("~ SRI*RIDRETH3_factor + SRI*RIAGENDR_factor + ", paste(cov_rest, collapse = " + "))), data = nd)
pred_one <- function(file) {
  r <- readRDS(file); cf <- r$coef; CF <- r$boot_coef
  stopifnot(setequal(colnames(X), names(cf))); Xm <- X[, names(cf)]
  if (anyNA(cf)) { bad <- names(cf)[is.na(cf)]; if (any(Xm[, bad] != 0)) stop("NA coefficient used by the reference person: ", paste(bad, collapse = ", ")); cf[is.na(cf)] <- 0 }
  CF[is.na(CF)] <- 0
  mu <- as.vector(Xm %*% cf); P <- Xm %*% t(CF)
  cbind(grid, input = r$input, calculator = r$calculator, n = r$n, bmi = exp(mu), lo = exp(apply(P, 1, quantile, .025)), hi = exp(apply(P, 1, quantile, .975)))
}
files <- list.files(m4dir, pattern = "^m4_.*[.]rds$", full.names = TRUE)
only <- Sys.getenv("ONLY", ""); if (nzchar(only)) files <- files[basename(files) %in% strsplit(only, ",")[[1]]]
pred <- rbindlist(lapply(files, pred_one)); fwrite(pred, file.path(m4dir, "m4_predictions.csv"))
cols <- c(Male = "#F8766D", Female = "#00BFC4"); shade <- c(Male = adjustcolor("#F8766D", 0.3), Female = adjustcolor("#00BFC4", 0.3))
labs <- c(published = "published SRI", raw = "released classification, unsmoothed", maj1h50 = "rolling majority, 1 h", ggir1h80 = "GGIR SRI1 smoothing, 1 h, 80%", bridge2 = "awakenings <30 min bridged", spt2 = "main sleep period only", `random sample` = "random sample")
calcs <- c(`published implementation` = "published implementation", SRI_hand = "Phillips formula", SRI_GGIR = "GGIR", SRI_sleepregE5 = "sleepreg, 5 days", SRI_sleepregE3 = "sleepreg, 3 days", PAXPREDM = "PAXPREDM", GGIR_SIB = "GGIR SIB", GGIR_SPT = "GGIR SPT", SWaN = "SWaN")
page <- function(p, ylim = c(20, 42)) {
  par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 3, 0))
  for (e in eth_lev) {
    plot(NA, xlim = c(0, 100), ylim = ylim, xlab = "SRI", ylab = "Predicted BMI (kg/m²)", main = if (e == "Other Race - Including Multi-Racial") "Other, incl. multiracial" else e, cex.main = 0.95, bty = "l", las = 1)
    for (s in sex_lev) { q <- p[eth == e & sex == s][order(SRI)]
      polygon(c(q$SRI, rev(q$SRI)), c(pmin(pmax(q$lo, ylim[1]), ylim[2]), rev(pmin(pmax(q$hi, ylim[1]), ylim[2]))), col = shade[[s]], border = NA)
      lines(q$SRI, q$bmi, col = cols[[s]], lwd = 2) }
    if (e == eth_lev[1]) legend("topright", legend = c(Male = "men", Female = "women")[sex_lev], col = cols[sex_lev], lwd = 2, bty = "n", cex = 0.9) }
  mtext(sprintf("Model 4: %s, %s (n = %d)", labs[[p$input[1]]], calcs[[p$calculator[1]]], p$n[1]), outer = TRUE, cex = 0.9, line = 1)
}
ord <- c("published", "raw", "maj1h50", "ggir1h80", "bridge2", "spt2", "random sample"); pred[, key := paste(input, calculator)]
keys <- unique(pred[order(match(input, ord), match(calculator, names(calcs))), key])
pdf(file.path(m4dir, "m4_predictions_all.pdf"), width = 11, height = 7); for (k in keys) page(pred[key == k]); invisible(dev.off())
for (k in c("published published implementation", "raw SRI_hand", "spt2 SRI_hand", "random sample PAXPREDM", "random sample SWaN")) if (k %in% keys) {
  png(file.path(m4dir, paste0("m4_pred_", gsub("[^A-Za-z0-9]", "_", k), ".png")), width = 11, height = 7, units = "in", res = 110); page(pred[key == k]); invisible(dev.off()) }
figdir <- Sys.getenv("FIG_DIR", ""); if (nzchar(figdir)) {                     # supplement figures: released series and main-sleep-period input, Phillips formula
  for (k in c("raw SRI_hand", "spt2 SRI_hand")) if (k %in% keys) { pdf(file.path(figdir, paste0("figS_model4_", sub(" SRI_hand", "", k), ".pdf")), width = 10, height = 6.8, pointsize = 13); page(pred[key == k]); invisible(dev.off()) } }
cat("pages:", length(keys), "| written to", m4dir, if (nzchar(figdir)) paste("and", figdir) else "", "\n")
