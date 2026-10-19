# Design-based (Taylor linearisation) 95% intervals of the fully adjusted and the sex-interaction model for comparison with the
# participant-bootstrap intervals of the paper (Supplementary Section 5). Env: SRI_PROJECT_ROOT (working-tree root).
suppressMessages({library(data.table); library(survey)})
base <- Sys.getenv("SRI_PROJECT_ROOT"); frag <- file.path(base, "01_recompute_SRI/intermediate/revision_2026-10_results/frag"); code <- file.path(base, "SRI_BMI_NHANES_2026_reanalysis/code/revision_2026-10")
rawdir <- file.path(base, "01_recompute_SRI/intermediate/revision_2026-10_results/raw_interim/r1000")
origdata <- file.path(base, "00_original_paper_materials", "data")
d  <- as.data.table(readRDS(file.path(origdata, "df_an1_imputed_MAIN.RDS")))
an <- as.data.table(readRDS(file.path(origdata, "Combined_data_sets_2_waves_analytic_data_set.RDS")))
d  <- merge(d, an[, .(SEQN, WTMEC2YR, SDMVPSU, SDMVSTRA)], by = "SEQN", all.x = TRUE)[, WTMEC4YR := WTMEC2YR / 2]
src <- readLines(file.path(code, "variants4.R")); s0 <- grep("^set.seed\\(3455\\)", src); s1 <- grep("^mf_table <- function", src) - 1; eval(parse(text = src[s0:s1]))
cat("strata:", uniqueN(d$SDMVSTRA), " PSUs:", uniqueN(paste(d$SDMVSTRA, d$SDMVPSU)), "\n")
taylor <- function(dd, svals, label) {
  dd <- copy(dd); set.seed(3455); dd[, SRI := winz(svals)]; dd <- dd[is.finite(SRI) & !is.na(BMXBMI) & !is.na(WTMEC4YR)]
  med <- dd[, median(SRI), by = .(q = cut(SRI, quantile(SRI, 0:5/5), include.lowest = TRUE, labels = 1:5))][order(q), V1]; dSRI <- med[5] - med[1]
  des <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, nest = TRUE, data = dd)
  f2 <- as.formula(paste0("log(BMXBMI) ~ SRI + ", paste(cov_full, collapse = " + "))); m <- svyglm(f2, design = des)
  b <- coef(m)["SRI"]; se <- sqrt(vcov(m)["SRI", "SRI"]); df <- degf(des); tq <- qt(.975, df)
  f3 <- as.formula(paste0("log(BMXBMI) ~ RIAGENDR_factor*SRI - SRI + SRI + ", paste(cov_full_nosex, collapse = " + "))); m3 <- svyglm(f3, design = des)
  cf <- coef(m3); V <- vcov(m3); ix <- grep("RIAGENDR_factor.*:SRI|SRI:RIAGENDR_factor", names(cf), value = TRUE)[1]
  bw <- cf["SRI"] + cf[ix]; sew <- sqrt(V["SRI", "SRI"] + V[ix, ix] + 2 * V["SRI", ix])
  cat(sprintf("%-34s n = %4d  M2 %.3f (%.3f, %.3f)   M3 women %.3f (%.3f, %.3f)   [design df = %d]\n", label, nrow(dd), exp(b * dSRI), exp((b - tq * se) * dSRI), exp((b + tq * se) * dSRI), exp(bw * dSRI), exp((bw - tq * sew) * dSRI), exp((bw + tq * sew) * dSRI), df)) }
cat("== Taylor-linearisation (design-based) 95% intervals; bootstrap in the paper: published 0.950 (0.930, 0.969) women 0.918 (0.891, 0.945); released 0.955 (0.936, 0.976) women 0.922 (0.902, 0.949)\n")
taylor(d, d$SRI_no_imputation_winsorized_rand, "published SRI")
g <- readRDS(file.path(frag, "na_raw.rds")); dm <- merge(d, g, by = "SEQN"); taylor(dm, dm$SRI_hand, "released series, Phillips formula")
g <- readRDS(file.path(frag, "na_spt2.rds")); dm <- merge(d, g, by = "SEQN"); taylor(dm, dm$SRI_hand, "main sleep period, Phillips formula")
source(file.path(code, "raw_sample.R")); res <- readRDS(file.path(rawdir, "raw_derivation_participants.rds")); ord <- fread(file.path(code, "analytic_random_order.tsv"))[, .(SEQN = seqn, rank)]
samp <- pick_sample(res, ord, "prefix", 1000, c(83452L, 69334L, 64226L))$samp; sw <- dcast(res[SEQN %in% samp], SEQN ~ deriv, value.var = "sri_A"); dr <- merge(d, sw, by = "SEQN")
dr <- dr[is.finite(PAXPREDM) & is.finite(GGIR_SIB) & is.finite(GGIR_SPT) & is.finite(SWaN)]
cat("== random sample; bootstrap in the paper: PAXPREDM 0.952 (0.907, 0.996) women 0.888 (0.819, 0.951); SWaN 0.993 (0.951, 1.025); GGIR SIB 1.004 (0.955, 1.041)\n")
for (s in c("PAXPREDM", "GGIR_SIB", "GGIR_SPT", "SWaN")) taylor(dr, dr[[s]], paste("random sample,", s))
