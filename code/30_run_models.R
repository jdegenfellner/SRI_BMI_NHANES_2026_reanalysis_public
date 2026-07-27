# =====================================================================
# PART 3: Re-estimate the published BMI-SRI association models for EACH
# SRI method (original/hand/GGIR/sleepreg) on the identical 7085-frame.
# Replicates model3 (fully adjusted, continuous), quintile MF Q5vsQ1,
# and sex/ethnicity interaction models. All survey-weighted (WTMEC4YR).
# =====================================================================
suppressMessages({library(data.table); library(survey)})
options(width = 170)

base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())
source(file.path(base, if (dir.exists(file.path(base, "code"))) "code" else "01_recompute_SRI", "00_paths.R"))

# ---- modeling frame (n=7085, all covariates + original SRI + BMI) ----
d  <- as.data.table(readRDS(file.path(origdata, "df_an1_imputed_MAIN.RDS")))
# survey weights from analytic dataset
an <- as.data.table(readRDS(file.path(origdata, "Combined_data_sets_2_waves_analytic_data_set.RDS")))
d  <- merge(d, an[, .(SEQN, WTMEC2YR, SDMVPSU, SDMVSTRA)], by = "SEQN", all.x = TRUE)
d[, WTMEC4YR := WTMEC2YR / 2]
# new SRI methods
sri <- as.data.table(readRDS(file.path(inter, "sri_all_methods.rds")))
d  <- merge(d, sri[, .(SEQN, SRI_hand, SRI_GGIR, SRI_sleepreg)], by = "SEQN", all.x = TRUE)
cat("modeling frame:", nrow(d), "rows\n")

# ---- winsorize-rand transform (same rule as paper: <-6.8 -> random-sleep noise)
set.seed(3455); SRIsim <- replicate(100000, -100 + 200*sum(rbinom(1440,1,.5))/1440)
m_sri <- mean(SRIsim); s_sri <- sd(SRIsim)
winz <- function(x){ i <- which(x < -6.8); x[i] <- rnorm(length(i), m_sri, s_sri); x }
d[, SRI_orig_w     := SRI_no_imputation_winsorized_rand]        # already winsorized in paper
d[, SRI_hand_w     := winz(SRI_hand)]
d[, SRI_GGIR_w     := winz(SRI_GGIR)]
d[, SRI_sleepreg_w := winz(SRI_sleepreg)]

methods <- c(orig="SRI_orig_w", hand="SRI_hand_w", GGIR="SRI_GGIR_w", sleepreg="SRI_sleepreg_w")

covars <- paste("RIDAGEYR","RIAGENDR_factor","RIDRETH3_factor","DMDEDUC2_factor",
                "INDHHIN2_factor","OCD150_factor","alcohol_consumption","smoking_status",
                "LBXVIDMS","DR1TKCAL_winsorized","Depression_score_PHQ_9","DMDMARTL_factor",
                "activity_level", sep=" + ")

run_one <- function(svar){
  dd <- d[!is.na(get(svar)) & !is.na(BMXBMI) & !is.na(WTMEC4YR)]
  dd[, SRIz := scale(get(svar))[,1]]                         # per-SD
  # quintiles (unweighted cutpoints) Q1..Q5
  qb <- quantile(dd[[svar]], probs = 0:5/5, na.rm = TRUE); qb[1] <- -Inf; qb[6] <- Inf
  dd[, SRIq := cut(get(svar), breaks = qb, labels = paste0("Q", 1:5), include.lowest = TRUE)]
  dd[, SRIq := relevel(factor(SRIq), ref = "Q1")]
  des <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC4YR, nest=TRUE, data=dd)

  # continuous, fully adjusted (= model3), per +1 SD; Wald z-test (survey df too small for t)
  m_cont <- svyglm(as.formula(paste0("log(BMXBMI) ~ SRIz + ", covars)), design=des)
  b <- unname(summary(m_cont)$coef["SRIz", ])
  cont <- c(MF_perSD = exp(b[1]), lo = exp(b[1]-1.96*b[2]), hi = exp(b[1]+1.96*b[2]),
            p = 2*pnorm(-abs(b[1]/b[2])))

  # quintile MF Q5 vs Q1, fully adjusted
  m_q <- svyglm(as.formula(paste0("log(BMXBMI) ~ SRIq + ", covars)), design=des)
  bq <- unname(summary(m_q)$coef["SRIqQ5", ])
  quint <- c(MF_Q5vQ1 = exp(bq[1]), lo = exp(bq[1]-1.96*bq[2]), hi = exp(bq[1]+1.96*bq[2]),
             p = 2*pnorm(-abs(bq[1]/bq[2])))

  # sex interaction (model4): SRIz * sex
  m_sex <- svyglm(as.formula(paste0("log(BMXBMI) ~ SRIz*RIAGENDR_factor + ",
            gsub("RIAGENDR_factor \\+ ", "", covars))), design=des)
  p_sex <- tryCatch(regTermTest(m_sex, ~SRIz:RIAGENDR_factor, df=Inf)$p[1], error=function(e) NA)

  # ethnicity interaction (model5): SRIz * ethnicity
  m_eth <- svyglm(as.formula(paste0("log(BMXBMI) ~ SRIz*RIDRETH3_factor + ",
            gsub("RIDRETH3_factor \\+ ", "", covars))), design=des)
  p_eth <- tryCatch(regTermTest(m_eth, ~SRIz:RIDRETH3_factor, df=Inf)$p[1], error=function(e) NA)

  list(n=nrow(dd), cont=cont, quint=quint, p_sex=p_sex, p_eth=p_eth,
       sri_mean=mean(dd[[svar]]), sri_sd=sd(dd[[svar]]))
}

cat("\n================ RESULTS BY SRI METHOD ================\n")
tab <- rbindlist(lapply(names(methods), function(nm){
  r <- run_one(methods[[nm]])
  data.table(method=nm, n=r$n,
    SRI_mean=round(r$sri_mean,2), SRI_sd=round(r$sri_sd,2),
    MF_perSD=round(r$cont["MF_perSD"],4),
    perSD_lo=round(r$cont["lo"],4), perSD_hi=round(r$cont["hi"],4),
    perSD_p=signif(r$cont["p"],3),
    MF_Q5vQ1=round(r$quint["MF_Q5vQ1"],4),
    Q5_lo=round(r$quint["lo"],4), Q5_hi=round(r$quint["hi"],4), Q5_p=signif(r$quint["p"],3),
    p_sexIA=signif(r$p_sex,3), p_ethIA=signif(r$p_eth,3))
}))
print(tab)
fwrite(tab, file.path(inter, "model_comparison.csv"))
saveRDS(tab, file.path(inter, "model_comparison.rds"))
cat("\nDONE part 3\n")
