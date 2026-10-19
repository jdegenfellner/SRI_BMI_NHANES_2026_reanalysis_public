# BMI models (published Table-3 recipe: survey-weighted log-BMI models, multiplicative factor Q1 -> Q5) for the SRI of each sleep/wake series
# derived from the 80 Hz raw signal, and for PAXPREDM on the same participants.
#
# Environment variables (all optional):
#   SRI_PROJECT_ROOT  working-tree root (contains 00_original_paper_materials/), default: working directory
#   IN_DIR            folder with raw_derivation_participants.rds (default: working directory)
#   OUT_DIR           output folder (default: IN_DIR)
#   ORDER_FILE, MODE, R_MAX, EXCLUDED   as in raw_derivation_summary.R (sample selection is shared through raw_sample.R)
#   SRI_COL           "sri_A" (common mask, main analysis; default) or "sri_B" (own missing minutes, sensitivity analysis; participants with >= 30% missing excluded)
#   NBOOT             bootstrap replications (default 200, as in the published recipe)
#
# The covariate sets, winsorisation rule, survey design and the effect size (difference of the median SRI between the fifth and the first quintile)
# are taken from variants4.R (same lines are evaluated). All series are fitted on the same participants and the same bootstrap samples, so the
# ratio of two multiplicative factors is formed within each replicate.
suppressMessages({library(data.table); library(survey)})
this <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE)   # Rscript encodes spaces in --file as ~+~
here <- if (is.na(this)) getwd() else dirname(normalizePath(this))
base   <- Sys.getenv("SRI_PROJECT_ROOT", getwd()); indir <- Sys.getenv("IN_DIR", getwd()); outdir <- Sys.getenv("OUT_DIR", indir)
ordf   <- Sys.getenv("ORDER_FILE", file.path(here, "analytic_random_order.tsv")); mode <- Sys.getenv("MODE", "prefix"); rmax <- as.numeric(Sys.getenv("R_MAX", "Inf"))
excl_f <- Sys.getenv("EXCLUDED", ""); scol <- Sys.getenv("SRI_COL", "sri_A"); nb <- as.integer(Sys.getenv("NBOOT", "200"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
source(file.path(here, "raw_sample.R"))

res <- readRDS(file.path(indir, "raw_derivation_participants.rds")); ord <- fread(ordf)[, .(SEQN = seqn, rank)]
excluded <- if (nzchar(excl_f) && file.exists(excl_f)) as.integer(readLines(excl_f)) else integer()
ps <- pick_sample(res, ord, mode, rmax, excluded); samp <- ps$samp
derivs <- c("PAXPREDM", "GGIR_SIB", "GGIR_SPT", "SWaN"); lab <- c(PAXPREDM = "PAXPREDM", GGIR_SIB = "GGIR SIB", GGIR_SPT = "GGIR SPT", SWaN = "SWaN")
mname <- c(PAXPREDM = "Paxpredm", GGIR_SIB = "GgirSib", GGIR_SPT = "GgirSpt", SWaN = "Swan")

# ---- data (as in variants4.R) ----------------------------------------------------------------------------------------------------
origdata <- file.path(base, "00_original_paper_materials", "data")
d  <- as.data.table(readRDS(file.path(origdata, "df_an1_imputed_MAIN.RDS")))
an <- as.data.table(readRDS(file.path(origdata, "Combined_data_sets_2_waves_analytic_data_set.RDS")))
d  <- merge(d, an[, .(SEQN, WTMEC2YR, SDMVPSU, SDMVSTRA)], by = "SEQN", all.x = TRUE)[, WTMEC4YR := WTMEC2YR / 2]
src <- readLines(file.path(here, "variants4.R"))
s0 <- grep("^set.seed\\(3455\\)", src); s1 <- grep("^mf_table <- function", src) - 1
eval(parse(text = src[s0:s1]))                        # SRIsim, winz, cov_full, beta_SRI, fit_set

sw <- dcast(res[SEQN %in% samp], SEQN ~ deriv, value.var = scol)
if (scol == "sri_B") { inc <- dcast(res[SEQN %in% samp], SEQN ~ deriv, value.var = "incl_B"); sw <- sw[SEQN %in% inc[GGIR_SIB & GGIR_SPT & SWaN, SEQN]] }
dd <- merge(d, sw, by = "SEQN")
for (s in derivs) { set.seed(20261017); dd[, (s) := winz(get(s))] }
dd <- dd[Reduce(`&`, lapply(derivs, function(s) is.finite(dd[[s]]))) & !is.na(BMXBMI) & !is.na(WTMEC4YR)]
n <- nrow(dd); cat(sprintf("sample %d, models on %d participants (SRI column %s)\n", length(samp), n, scol))

eff <- lapply(derivs, function(s) { x <- dd[[s]]; med <- tapply(x, cut(x, quantile(x, 0:5 / 5), include.lowest = TRUE, labels = 1:5), median); list(d = unname(med[5] - med[1]), q1 = unname(med[1]), q5 = unname(med[5])) }); names(eff) <- derivs
fit_series <- function(data, s) { data <- copy(data); data[, SRI := get(s)]; exp(fit_set(data) * eff[[s]]$d) }
pnames <- c("MA", "MB", "MCref", "MCother")       # age+sex model, fully adjusted model, fully adjusted model by sex (reference group, other group)
pt <- sapply(derivs, function(s) fit_series(dd, s)); rownames(pt) <- pnames
set.seed(11)
B <- array(NA_real_, c(nb, length(derivs), 4), dimnames = list(NULL, derivs, pnames))
for (b in seq_len(nb)) { bs <- dd[sample(n, replace = TRUE)]
  for (s in derivs) B[b, s, ] <- tryCatch(fit_series(bs, s), error = function(e) rep(NA_real_, 4))
  if (b %% 25 == 0) cat(sprintf("  bootstrap %d/%d\n", b, nb)) }
q <- function(x) quantile(x, c(.025, .975), na.rm = TRUE)
out <- rbindlist(lapply(derivs, function(s) rbindlist(lapply(pnames, function(p) { ci <- q(B[, s, p])
  data.table(series = lab[[s]], model = p, n = n, dSRI = eff[[s]]$d, medQ1 = eff[[s]]$q1, medQ5 = eff[[s]]$q5, mf = pt[p, s], lo = ci[1], hi = ci[2]) }))))
rat <- rbindlist(lapply(derivs[-1], function(s) rbindlist(lapply(pnames, function(p) { r <- B[, s, p] / B[, "PAXPREDM", p]; ci <- q(r)
  data.table(series = lab[[s]], model = p, ratio_to_PAXPREDM = pt[p, s] / pt[p, "PAXPREDM"], lo = ci[1], hi = ci[2]) }))))
fwrite(out, file.path(outdir, paste0("raw_models_", scol, ".csv"))); fwrite(rat, file.path(outdir, paste0("raw_models_ratio_", scol, ".csv")))
saveRDS(list(out = out, ratio = rat, boot = B, n = n), file.path(outdir, paste0("raw_models_", scol, ".rds")))

# ---- LaTeX macros -----------------------------------------------------------------------------------------------------------------
suffix <- if (scol == "sri_B") "B" else ""
mac <- character(); add <- function(nm, v) mac <<- c(mac, sprintf("\\newcommand{\\%s}{%s}", nm, v))
fsgn <- function(x) sub("^-0\\.0$", "0.0", sprintf("%+.1f", x))   # signed one-decimal percent, without a negative zero
add(paste0("rawModelN", suffix), n)
for (s in derivs) for (p in pnames) { r <- out[series == lab[[s]] & model == p]; nm <- paste0("raw", mname[[s]], p, suffix)
  add(paste0(nm, "Est"), sprintf("%.3f", r$mf)); add(paste0(nm, "Lo"), sprintf("%.3f", r$lo)); add(paste0(nm, "Hi"), sprintf("%.3f", r$hi)); add(paste0(nm, "DSri"), sprintf("%.1f", r$dSRI))
  pt10 <- 100 * (exp(log(r$mf) / r$dSRI * 10) - 1); ci10 <- 100 * (exp(log(q(B[, s, p])) / r$dSRI * 10) - 1)   # percent change in BMI per 10 SRI points
  add(paste0(nm, "PerTenEst"), fsgn(pt10)); add(paste0(nm, "PerTenLo"), fsgn(ci10[1])); add(paste0(nm, "PerTenHi"), fsgn(ci10[2])) }
for (s in derivs[-1]) for (p in pnames) { r <- rat[series == lab[[s]] & model == p]; nm <- paste0("raw", mname[[s]], p, "Ratio", suffix)
  add(paste0(nm, "Est"), sprintf("%.3f", r$ratio_to_PAXPREDM)); add(paste0(nm, "Lo"), sprintf("%.3f", r$lo)); add(paste0(nm, "Hi"), sprintf("%.3f", r$hi)) }
writeLines(c("% generated by raw_derivation_models.R; do not edit", mac), file.path(outdir, paste0("raw_models_numbers", suffix, ".tex")))
cat("\n== multiplicative factors Q1 -> Q5 (95% bootstrap interval) ==\n")
print(out[, .(series, model, n, dSRI = round(dSRI, 1), MF = sprintf("%.3f (%.3f, %.3f)", mf, lo, hi))])
cat("\n== ratio of the multiplicative factor to PAXPREDM (same participants, same bootstrap samples) ==\n")
print(rat[, .(series, model, ratio = sprintf("%.3f (%.3f, %.3f)", ratio_to_PAXPREDM, lo, hi))])
