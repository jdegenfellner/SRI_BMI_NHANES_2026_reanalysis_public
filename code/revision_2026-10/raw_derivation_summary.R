# Summary of the sleep/wake series derived from the 80 Hz raw signal (GGIR SIB, GGIR SPT, SWaN) against the released classification (PAXPREDM).
# Input: raw_derivation_participants.rds and raw_derivation_bouts.rds (from raw_derivation_participants.R).
#
# Environment variables (all optional):
#   IN_DIR       folder with the two RDS files (default: working directory)
#   OUT_DIR      output folder (default: IN_DIR)
#   ORDER_FILE   analytic_random_order.tsv (columns rank, seqn, cycle); default: analytic_random_order.tsv next to this script
#   MODE         "prefix" (default): sample = participants with rank <= R whose derived series are complete, R = largest rank such that
#                every participant with a smaller or equal rank is complete or listed in EXCLUDED (random sample of the analytic sample);
#                "all": every participant with all four series (final analysis)
#   R_MAX        upper limit for R in prefix mode (default: no limit)
#   STEPS        comma-separated prefix sizes for the stability table (default 500,1000,1500,...,7000,7085)
#   EXCLUDED     optional text file with SEQN that failed in the derivation and are reported as missing (one per line)
#   NBOOT        bootstrap replications (default 2000)
#
# Estimates are unweighted comparisons of the series. Intervals are percentile intervals of a participant-level bootstrap
# (participants are resampled; every statistic is recomputed in each replicate, differences are formed within participant).
suppressMessages(library(data.table))
this <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE)   # Rscript encodes spaces in --file as ~+~
here <- if (is.na(this)) getwd() else dirname(normalizePath(this))
indir  <- Sys.getenv("IN_DIR", getwd()); outdir <- Sys.getenv("OUT_DIR", indir)
ordf   <- Sys.getenv("ORDER_FILE", file.path(here, "analytic_random_order.tsv"))
mode   <- Sys.getenv("MODE", "prefix"); rmax_env <- as.numeric(Sys.getenv("R_MAX", "Inf"))
steps  <- as.integer(strsplit(Sys.getenv("STEPS", paste(c(seq(500, 7000, 500), 7085), collapse = ",")), ",")[[1]])
nboot  <- as.integer(Sys.getenv("NBOOT", "2000")); excl_f <- Sys.getenv("EXCLUDED", "")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

res <- readRDS(file.path(indir, "raw_derivation_participants.rds")); bouts <- readRDS(file.path(indir, "raw_derivation_bouts.rds"))
ord <- fread(ordf)[, .(SEQN = seqn, rank)]
derivs <- c("PAXPREDM", "GGIR_SIB", "GGIR_SPT", "SWaN"); der <- derivs[-1]
lab <- c(PAXPREDM = "PAXPREDM", GGIR_SIB = "GGIR SIB", GGIR_SPT = "GGIR SPT", SWaN = "SWaN")
mname <- c(PAXPREDM = "Paxpredm", GGIR_SIB = "GgirSib", GGIR_SPT = "GgirSpt", SWaN = "Swan")

# ---- complete participants and sample --------------------------------------------------------------------------------------
source(file.path(here, "raw_sample.R"))
excluded <- if (nzchar(excl_f) && file.exists(excl_f)) as.integer(readLines(excl_f)) else integer()
ps <- pick_sample(res, ord, mode, rmax_env, excluded); samp <- ps$samp; complete <- ps$complete; ord <- ps$ord; R <- ps$R
cat(sprintf("mode %s: %d participants with all four series; sample %d%s\n", mode, length(complete), length(samp), if (mode == "prefix") sprintf(" (prefix R = %d, %d excluded as failed)", R, ps$n_failed) else ""))
if (length(samp) < 10) stop("fewer than 10 participants in the sample")

# ---- statistics ------------------------------------------------------------------------------------------------------------------
kappa_of <- function(n11, n10, n01, n00) { N <- n11 + n10 + n01 + n00; po <- (n11 + n00) / N
  pe <- ((n11 + n10) * (n11 + n01) + (n01 + n00) * (n10 + n00)) / N^2; (po - pe) / (1 - pe) }
# per-participant wide matrices (rows = participants in `ids`)
wide <- function(ids, col) { x <- dcast(res[SEQN %in% ids], SEQN ~ deriv, value.var = col); x[match(ids, SEQN)] }
point_and_boot <- function(ids, B = nboot, seed = 20261017) {
  n <- length(ids); S <- wide(ids, "sri_A"); T <- wide(ids, "trans_day_A"); SL <- wide(ids, "sleep_min_day_A")
  cn <- lapply(c("all", "ism"), function(p) lapply(der, function(d) {
    x <- res[SEQN %in% ids & deriv == d]; x <- x[match(ids, SEQN)]
    as.matrix(x[, paste0(p, c("_N", "_agree", "_dmiss", "_n11", "_n10", "_n01", "_n00")), with = FALSE]) })); names(cn) <- c("all", "ism"); cn$non <- Map(function(a, b) a - b, cn$all, cn$ism)
  stat <- function(w) {
    out <- list()
    for (d in derivs) { out[[paste0("sri_", d)]] <- sum(w * S[[d]]) / sum(w); out[[paste0("trans_", d)]] <- sum(w * T[[d]]) / sum(w); out[[paste0("sleep_", d)]] <- sum(w * SL[[d]]) / sum(w) }
    for (d in der) out[[paste0("diff_", d)]] <- sum(w * (S[[d]] - S$PAXPREDM)) / sum(w)
    for (p in c("all", "ism", "non")) for (k in seq_along(der)) { m <- colSums(w * cn[[p]][[k]]); d <- der[k]
      out[[paste0("agree_", p, "_", d)]] <- m[2] / m[1]
      out[[paste0("kappa_", p, "_", d)]] <- kappa_of(m[4], m[5], m[6], m[7]) }
    unlist(lapply(out, unname)) }
  est <- stat(rep(1, n)); set.seed(seed)
  bs <- vapply(seq_len(B), function(i) stat(tabulate(sample.int(n, n, TRUE), n)), numeric(length(est)))
  data.table(stat = names(est), est = est, lo = apply(bs, 1, quantile, 0.025), hi = apply(bs, 1, quantile, 0.975))
}
main <- point_and_boot(samp)
gs <- function(nm, tab = main) tab[stat == nm]

# ---- tables ----------------------------------------------------------------------------------------------------------------------
sri_tab <- rbindlist(lapply(derivs, function(d) { x <- res[SEQN %in% samp & deriv == d]
  data.table(series = lab[[d]], n = nrow(x), mean = main[stat == paste0("sri_", d), est], sd = sd(x$sri_A), median = median(x$sri_A), q1 = quantile(x$sri_A, .25), q3 = quantile(x$sri_A, .75),
             diff_mean = if (d == "PAXPREDM") NA else gs(paste0("diff_", d))$est, diff_lo = if (d == "PAXPREDM") NA else gs(paste0("diff_", d))$lo, diff_hi = if (d == "PAXPREDM") NA else gs(paste0("diff_", d))$hi,
             trans_day = main[stat == paste0("trans_", d), est], sleep_min_day = main[stat == paste0("sleep_", d), est],
             miss_in_mask = mean(x$miss_in_mask), sri_B_mean = mean(x[incl_B == TRUE, sri_B]), n_incl_B = sum(x$incl_B), pct_excl_B = 100 * mean(!x$incl_B)) }))
agree_tab <- rbindlist(lapply(der, function(d) { x <- res[SEQN %in% samp & deriv == d]
  rbindlist(lapply(c("all", "ism", "non"), function(p) { N <- if (p == "non") sum(x$all_N) - sum(x$ism_N) else sum(x[[paste0(p, "_N")]])
    data.table(series = lab[[d]], minutes = p, N = N, pct_agree = 100 * main[stat == paste0("agree_", p, "_", d), est], pct_derived_missing = 100 * (if (p == "non") sum(x$all_dmiss) - sum(x$ism_dmiss) else sum(x[[paste0(p, "_dmiss")]])) / N,
               kappa = main[stat == paste0("kappa_", p, "_", d), est], kappa_lo = main[stat == paste0("kappa_", p, "_", d), lo], kappa_hi = main[stat == paste0("kappa_", p, "_", d), hi])
  })) }))
b <- bouts[SEQN %in% samp]
bout_tab <- b[, .(n_runs = .N, runs_per_participant = .N / length(samp), median_min = as.numeric(median(len)), q1 = as.numeric(quantile(len, .25)), q3 = as.numeric(quantile(len, .75)), mean_min = mean(len)), by = .(series = lab[deriv], state = fifelse(state == 1, "sleep", "wake"))]
setorder(bout_tab, series, state)
fwrite(sri_tab, file.path(outdir, "raw_summary_sri.csv")); fwrite(agree_tab, file.path(outdir, "raw_summary_agreement.csv")); fwrite(bout_tab, file.path(outdir, "raw_summary_bouts.csv"))
fwrite(main, file.path(outdir, "raw_summary_estimates.csv"))

# ---- stability over prefix sizes -----------------------------------------------------------------------------------------------
stab <- NULL
if (mode == "prefix") {
  st <- sort(unique(c(steps[steps <= R], R))); st <- st[st >= as.integer(Sys.getenv("MINSTEP", "100"))]
  stab <- rbindlist(lapply(st, function(r) { ids <- ord[rank <= r & complete, SEQN]; x <- point_and_boot(ids, B = min(nboot, 1000)); x[, `:=`(r = r, n = length(ids))]; x }))
  fwrite(stab, file.path(outdir, "raw_summary_stability.csv"))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    library(ggplot2)
    sel <- stab[stat %in% c(paste0("diff_", der), paste0("kappa_all_", der), paste0("sri_", derivs))]
    sel[, panel := fifelse(grepl("^diff_", stat), "Mean SRI difference to PAXPREDM", fifelse(grepl("^kappa_", stat), "Pooled kappa with PAXPREDM", "Mean SRI"))]
    sel[, series := lab[sub("^(diff|kappa_all|sri)_", "", stat)]]
    g <- ggplot(sel, aes(r, est, colour = series, fill = series)) + geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line() + geom_point(size = 1) +
      facet_wrap(~panel, scales = "free_y", ncol = 1) + labs(x = "Participants in the random sample (prefix of the fixed random order)", y = NULL, colour = NULL, fill = NULL) + theme_minimal(base_size = 10) + theme(legend.position = "bottom")
    ggsave(file.path(outdir, "raw_summary_stability.pdf"), g, width = 6.5, height = 8)
  }
}

# ---- LaTeX macros -----------------------------------------------------------------------------------------------------------------
f1 <- function(x) sprintf("%.1f", x); f2 <- function(x) sprintf("%.2f", x); f3 <- function(x) sprintf("%.3f", x)
mac <- character()
add <- function(name, value) mac <<- c(mac, sprintf("\\newcommand{\\%s}{%s}", name, value))
add("rawMode", mode); add("rawNsample", length(samp)); add("rawNcomplete", length(complete)); if (mode == "prefix") add("rawPrefixR", R)
for (d in derivs) { s <- sri_tab[series == lab[[d]]]; nm <- paste0("raw", mname[[d]])
  add(paste0(nm, "SriMean"), f1(s$mean)); add(paste0(nm, "SriSd"), f1(s$sd)); add(paste0(nm, "SriMedian"), f1(s$median)); add(paste0(nm, "TransDay"), f1(s$trans_day)); add(paste0(nm, "SleepMinDay"), sprintf("%.0f", s$sleep_min_day))
  add(paste0(nm, "MissMask"), f1(100 * s$miss_in_mask)); add(paste0(nm, "SriBMean"), f1(s$sri_B_mean)); add(paste0(nm, "PctExclB"), f1(s$pct_excl_B))
  if (d != "PAXPREDM") { add(paste0(nm, "DiffMean"), f1(s$diff_mean)); add(paste0(nm, "DiffLo"), f1(s$diff_lo)); add(paste0(nm, "DiffHi"), f1(s$diff_hi))
    for (p in c("all", "ism", "non")) { a <- agree_tab[series == lab[[d]] & minutes == p]; pn <- c(all = "All", ism = "Ism", non = "Non")[[p]]
      add(paste0(nm, "Agree", pn), f1(a$pct_agree)); add(paste0(nm, "Kappa", pn), f2(a$kappa)); add(paste0(nm, "Kappa", pn, "Lo"), f2(a$kappa_lo)); add(paste0(nm, "Kappa", pn, "Hi"), f2(a$kappa_hi)); add(paste0(nm, "Dmiss", pn), f1(a$pct_derived_missing)) } } }
for (i in seq_len(nrow(bout_tab))) { x <- bout_tab[i]; nm <- paste0("raw", mname[[names(lab)[match(x$series, lab)]]], "Bout", if (x$state == "sleep") "Sleep" else "Wake")
  add(paste0(nm, "Median"), f1(x$median_min)); add(paste0(nm, "PerPerson"), f1(x$runs_per_participant)) }
writeLines(c("% generated by raw_derivation_summary.R; do not edit", mac), file.path(outdir, "raw_derivation_numbers.tex"))
cat("\n== SRI ==\n"); print(sri_tab[, .(series, n, mean = round(mean, 1), sd = round(sd, 1), diff = round(diff_mean, 1), lo = round(diff_lo, 1), hi = round(diff_hi, 1), trans = round(trans_day, 1), miss = round(miss_in_mask, 3), sriB = round(sri_B_mean, 1), exclB = round(pct_excl_B, 1))])
cat("\n== agreement ==\n"); print(agree_tab[, .(series, minutes, N, agree = round(pct_agree, 1), dmiss = round(pct_derived_missing, 1), kappa = round(kappa, 2), lo = round(kappa_lo, 2), hi = round(kappa_hi, 2))])
cat("\n== bouts (minutes) ==\n"); print(bout_tab[, .(series, state, runs_pp = round(runs_per_participant, 1), median = median_min, q1, q3)])
if (!is.null(stab)) { cat("\n== stability (mean difference to PAXPREDM) ==\n"); print(dcast(stab[grepl("^diff_", stat)], r + n ~ stat, value.var = "est")[, lapply(.SD, function(x) round(x, 2))]) }
cat("macros:", length(mac), "->", file.path(outdir, "raw_derivation_numbers.tex"), "\n")
