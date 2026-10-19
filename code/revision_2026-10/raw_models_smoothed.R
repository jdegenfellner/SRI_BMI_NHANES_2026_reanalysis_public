# Models 1 to 3 (published recipe) for the four smoothed versions of the released classification on the same participants and the same
# bootstrap resamples as raw_derivation_models.R, so that Table 4 of the main text shows the released classification, its smoothed
# versions and the series derived from the raw signal side by side. The smoothed SRIs are the Phillips-formula values of the input
# variants (frag/na_<variant>.rds, column SRI_hand, all minutes of the rebuilt series).
# Environment variables as in raw_derivation_models.R (SRI_PROJECT_ROOT, IN_DIR, OUT_DIR, ORDER_FILE, MODE, R_MAX, EXCLUDED, NBOOT) plus
# FRAG_DIR (folder with na_<variant>.rds; default <SRI_PROJECT_ROOT>/01_recompute_SRI/intermediate/revision_2026-10_results/frag).
suppressMessages({library(data.table); library(survey)})
this <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE)
here <- if (is.na(this)) getwd() else dirname(normalizePath(this))
base   <- Sys.getenv("SRI_PROJECT_ROOT", getwd()); indir <- Sys.getenv("IN_DIR", getwd()); outdir <- Sys.getenv("OUT_DIR", indir)
fragd  <- Sys.getenv("FRAG_DIR", file.path(base, "01_recompute_SRI/intermediate/revision_2026-10_results/frag"))
ordf   <- Sys.getenv("ORDER_FILE", file.path(here, "analytic_random_order.tsv")); mode <- Sys.getenv("MODE", "prefix"); rmax <- as.numeric(Sys.getenv("R_MAX", "Inf"))
excl_f <- Sys.getenv("EXCLUDED", ""); nb <- as.integer(Sys.getenv("NBOOT", "200"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE); source(file.path(here, "raw_sample.R"))
res <- readRDS(file.path(indir, "raw_derivation_participants.rds")); ord <- fread(ordf)[, .(SEQN = seqn, rank)]
excluded <- if (nzchar(excl_f) && file.exists(excl_f)) as.integer(readLines(excl_f)) else as.integer(strsplit(Sys.getenv("EXCLUDED_SEQN", ""), ",")[[1]])
ps <- pick_sample(res, ord, mode, rmax, excluded); samp <- ps$samp
derivs <- c("PAXPREDM", "GGIR_SIB", "GGIR_SPT", "SWaN")
variants <- c(maj1h50 = "Maj", ggir1h80 = "Ggirsm", bridge2 = "Bridge", spt2 = "Spt"); lab <- c(maj1h50 = "Rolling majority", ggir1h80 = "GGIR SRI1 smoothing", bridge2 = "Bridged", spt2 = "Main sleep period")
origdata <- file.path(base, "00_original_paper_materials", "data")
d  <- as.data.table(readRDS(file.path(origdata, "df_an1_imputed_MAIN.RDS")))
an <- as.data.table(readRDS(file.path(origdata, "Combined_data_sets_2_waves_analytic_data_set.RDS")))
d  <- merge(d, an[, .(SEQN, WTMEC2YR, SDMVPSU, SDMVSTRA)], by = "SEQN", all.x = TRUE)[, WTMEC4YR := WTMEC2YR / 2]
src <- readLines(file.path(here, "variants4.R")); s0 <- grep("^set.seed\\(3455\\)", src); s1 <- grep("^mf_table <- function", src) - 1; eval(parse(text = src[s0:s1]))
sw <- dcast(res[SEQN %in% samp], SEQN ~ deriv, value.var = "sri_A"); dd <- merge(d, sw, by = "SEQN")       # same construction as raw_derivation_models.R
for (s in derivs) { set.seed(20261017); dd[, (s) := winz(get(s))] }
dd <- dd[Reduce(`&`, lapply(derivs, function(s) is.finite(dd[[s]]))) & !is.na(BMXBMI) & !is.na(WTMEC4YR)]
for (v in names(variants)) { g <- readRDS(file.path(fragd, paste0("na_", v, ".rds")))[, .(SEQN, SRI_hand)]; dd[, (v) := g[dd, on = "SEQN", SRI_hand]]; set.seed(20261017); dd[, (v) := winz(get(v))] }
stopifnot(all(sapply(names(variants), function(v) all(is.finite(dd[[v]])))))
n <- nrow(dd); cat(sprintf("sample %d, models on %d participants\n", length(samp), n))
series <- c("PAXPREDM", names(variants))
eff <- lapply(series, function(s) { x <- dd[[s]]; med <- tapply(x, cut(x, quantile(x, 0:5 / 5), include.lowest = TRUE, labels = 1:5), median); list(d = unname(med[5] - med[1])) }); names(eff) <- series
fit_series <- function(data, s) { data <- copy(data); data[, SRI := get(s)]; exp(fit_set(data) * eff[[s]]$d) }
pnames <- c("MA", "MB", "MCref", "MCother"); pt <- sapply(series, function(s) fit_series(dd, s)); rownames(pt) <- pnames
set.seed(11); B <- array(NA_real_, c(nb, length(series), 4), dimnames = list(NULL, series, pnames))
for (b in seq_len(nb)) { bs <- dd[sample(n, replace = TRUE)]                                              # identical resamples to raw_derivation_models.R
  for (s in series) B[b, s, ] <- tryCatch(fit_series(bs, s), error = function(e) rep(NA_real_, 4))
  if (b %% 50 == 0) cat(sprintf("  bootstrap %d/%d\n", b, nb)) }
q <- function(x) quantile(x, c(.025, .975), na.rm = TRUE)
out <- rbindlist(lapply(series, function(s) rbindlist(lapply(pnames, function(p) { ci <- q(B[, s, p]); data.table(series = s, model = p, n = n, dSRI = eff[[s]]$d, mf = pt[p, s], lo = ci[1], hi = ci[2]) }))))
rat <- rbindlist(lapply(names(variants), function(s) rbindlist(lapply(pnames, function(p) { r <- B[, s, p] / B[, "PAXPREDM", p]; ci <- q(r); data.table(series = s, model = p, ratio = pt[p, s] / pt[p, "PAXPREDM"], lo = ci[1], hi = ci[2]) }))))
fwrite(out, file.path(outdir, "raw_models_smoothed.csv")); fwrite(rat, file.path(outdir, "raw_models_smoothed_ratio.csv")); saveRDS(list(out = out, ratio = rat, boot = B, n = n), file.path(outdir, "raw_models_smoothed.rds"))
mac <- character(); add <- function(nm, v) mac <<- c(mac, sprintf("\\newcommand{\\%s}{%s}", nm, v))
for (v in names(variants)) for (p in pnames) { r <- out[series == v & model == p]; nm <- paste0("raw", variants[[v]], p)
  add(paste0(nm, "Est"), sprintf("%.3f", r$mf)); add(paste0(nm, "Lo"), sprintf("%.3f", r$lo)); add(paste0(nm, "Hi"), sprintf("%.3f", r$hi)); add(paste0(nm, "DSri"), sprintf("%.1f", r$dSRI))
  rr <- rat[series == v & model == p]; add(paste0(nm, "RatioEst"), sprintf("%.3f", rr$ratio)); add(paste0(nm, "RatioLo"), sprintf("%.3f", rr$lo)); add(paste0(nm, "RatioHi"), sprintf("%.3f", rr$hi)) }
mb <- out[model == "MB" & series != "PAXPREDM"]; add("rawSmoothedMBMin", sprintf("%.3f", min(mb$mf))); add("rawSmoothedMBMax", sprintf("%.3f", max(mb$mf)))
writeLines(c("% generated by raw_models_smoothed.R; do not edit", mac), file.path(outdir, "raw_models_smoothed_numbers.tex"))
cat("\n== check: PAXPREDM factors (must equal raw_models_numbers.tex) ==\n"); print(out[series == "PAXPREDM", .(model, MF = sprintf("%.3f (%.3f, %.3f)", mf, lo, hi))])
cat("\n== smoothed versions, Model 2 ==\n"); print(out[model == "MB", .(series, dSRI = round(dSRI, 1), MF = sprintf("%.3f (%.3f, %.3f)", mf, lo, hi))]); print(rat[model == "MB", .(series, ratio = sprintf("%.3f (%.3f, %.3f)", ratio, lo, hi))])
