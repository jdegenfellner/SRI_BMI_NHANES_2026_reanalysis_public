# Sensitivity check: the earlier study modelled log(BMXBMI_winsorized) (BMI winsorised at the 1st and 99th percentiles) in its
# Models 2 to 4; the re-analysis models log(BMXBMI). This script refits Models 1 to 4 with both outcomes on the published SRI and on
# the Phillips formula of the rebuilt released series and prints the differences of the multiplication factors (point estimates).
# Env: SRI_PROJECT_ROOT (working-tree root), FRAG_DIR (folder with na_raw.rds), OUT_DIR (default: working directory).
suppressMessages({library(data.table); library(survey)})
this <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE)
here <- if (is.na(this)) getwd() else dirname(normalizePath(this))
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd()); outdir <- Sys.getenv("OUT_DIR", getwd())
fragdir <- Sys.getenv("FRAG_DIR", file.path(base, "01_recompute_SRI/intermediate/revision_2026-10_results/frag"))
origdata <- file.path(base, "00_original_paper_materials", "data")
d  <- as.data.table(readRDS(file.path(origdata, "df_an1_imputed_MAIN.RDS")))
an <- as.data.table(readRDS(file.path(origdata, "Combined_data_sets_2_waves_analytic_data_set.RDS")))
d  <- merge(d, an[, .(SEQN, WTMEC2YR, SDMVPSU, SDMVSTRA)], by = "SEQN", all.x = TRUE)[, WTMEC4YR := WTMEC2YR / 2]
src <- readLines(file.path(here, "variants4.R")); s0 <- grep("^set.seed\\(3455\\)", src); s1 <- grep("^mf_table <- function", src) - 1
eval(parse(text = src[s0:s1]))                        # SRIsim, winz, cov_full, cov_full_nosex, beta_SRI, fit_set
d <- merge(d, readRDS(file.path(fragdir, "na_raw.rds"))[, .(SEQN, SRI_hand)], by = "SEQN", all.x = TRUE)
cov_rest <- setdiff(cov_full, c("RIAGENDR_factor", "RIDRETH3_factor"))
prep <- function(dd, svals) { dd <- copy(dd); set.seed(1); dd[, SRI := winz(svals)]; dd <- dd[is.finite(SRI) & !is.na(BMXBMI) & !is.na(WTMEC4YR)]
  med <- dd[, median(SRI), by = .(q = cut(SRI, quantile(SRI, 0:5 / 5), include.lowest = TRUE, labels = 1:5))][order(q), V1]
  list(dd = dd, dS = med[5] - med[1], des = svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, nest = TRUE, data = dd)) }
m123 <- function(p, outcome) { bf <- function(rhs) coef(svyglm(as.formula(paste0("log(", outcome, ") ~ SRI + ", rhs)), design = p$des))
  b1 <- bf("RIDAGEYR + RIAGENDR_factor")["SRI"]; b2 <- bf(paste(cov_full, collapse = " + "))["SRI"]
  c3 <- bf(paste0("RIAGENDR_factor*SRI - SRI + SRI + ", paste(cov_full_nosex, collapse = " + "))); ix <- grep("RIAGENDR_factor.*:SRI|SRI:RIAGENDR_factor", names(c3), value = TRUE)
  exp(c(M1 = unname(b1), M2 = unname(b2), M3_men = unname(c3["SRI"]), M3_women = unname(c3["SRI"] + c3[ix[1]])) * p$dS) }
m4 <- function(p, outcome) { co <- coef(svyglm(as.formula(paste0("log(", outcome, ") ~ SRI*RIDRETH3_factor + SRI*RIAGENDR_factor + ", paste(cov_rest, collapse = " + "))), design = p$des))
  out <- c(); for (e in levels(p$dd$RIDRETH3_factor)) for (s in levels(p$dd$RIAGENDR_factor)) { b <- co["SRI"]
    ie <- grep(paste0("^SRI:RIDRETH3_factor", e, "$|^RIDRETH3_factor", e, ":SRI$"), names(co)); is <- grep(paste0("^SRI:RIAGENDR_factor", s, "$|^RIAGENDR_factor", s, ":SRI$"), names(co))
    if (length(ie)) b <- b + co[ie]; if (length(is)) b <- b + co[is]; out[paste(e, s, sep = " / ")] <- exp(b * p$dS) }; out }
rows <- list()
for (s in c("SRI_no_imputation", "SRI_hand")) { p <- prep(d, d[[s]])
  a <- c(m123(p, "BMXBMI"), m4(p, "BMXBMI")); b <- c(m123(p, "BMXBMI_winsorized"), m4(p, "BMXBMI_winsorized"))
  rows[[s]] <- data.table(sri = s, model = names(a), logBMI = round(a, 4), logBMI_winsorised = round(b, 4), diff = round(b - a, 4))
  cat("\n==", s, ": participants with SRI < -6.8 replaced:", sum(d[[s]] < -6.8, na.rm = TRUE), "\n"); print(rows[[s]])
  cat("max |diff| Models 1-3:", round(max(abs(b[1:4] - a[1:4])), 4), "  Model 4:", round(max(abs(b[-(1:4)] - a[-(1:4)])), 4), "\n") }
fwrite(rbindlist(rows), file.path(outdir, "check_bmi_winsorisation.csv"))
