# =====================================================================
# Model-assumption checks for the NEW models (the fully adjusted Model 2 of
# Degenfellner et al. 2026, Table 3, re-fitted with each SRI calculator):
#   linearity of the SRI term, homoscedasticity, normality of residuals,
#   multicollinearity. Quantitative table + diagnostic plots, all four
#   calculators. (Unweighted lm, matching the original assumption checks.)
# =====================================================================
suppressMessages({library(data.table); library(ggplot2); library(patchwork)
                  library(lmtest); library(car); library(moments)})
options(width=170)
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())
source(file.path(base, if (dir.exists(file.path(base, "code"))) "code" else "01_recompute_SRI", "00_paths.R"))
# sdir provided by 00_paths.R

d  <- as.data.table(readRDS(file.path(origdata, "df_an1_imputed_MAIN.RDS")))
sri<- as.data.table(readRDS(file.path(inter,"sri_all_methods.rds")))[,.(SEQN,SRI_GGIR,SRI_sleepreg)]
d  <- merge(d, sri, by="SEQN", all.x=TRUE)
covs <- "RIDAGEYR + RIAGENDR_factor + RIDRETH3_factor + DMDEDUC2_factor + INDHHIN2_factor + OCD150_factor + alcohol_consumption + smoking_status + LBXVIDMS + DR1TKCAL_winsorized + Depression_score_PHQ_9 + DMDMARTL_factor + activity_level"
calcs <- c(Published="SRI_no_imputation_winsorized_rand", GGIR="SRI_GGIR", sleepreg="SRI_sleepreg")

rows<-list(); diaglist<-list()
for(nm in names(calcs)){ v<-calcs[[nm]]
  dd<-d[is.finite(get(v))&is.finite(BMXBMI)]
  m <-lm(as.formula(paste0("log(BMXBMI) ~ ",v," + ",covs)), data=dd)
  r <-residuals(m)
  # linearity of the SRI term: add quadratic, test improvement
  m2<-lm(as.formula(paste0("log(BMXBMI) ~ poly(",v,",2) + ",covs)), data=dd)
  p_nonlin<-anova(m,m2)$`Pr(>F)`[2]
  bp<-tryCatch(bptest(m)$p.value, error=function(e)NA)             # homoscedasticity
  vmax<-tryCatch(max(car::vif(m)[,1]), error=function(e) max(car::vif(m)))  # multicollinearity (GVIF or VIF)
  rows[[nm]]<-data.table(calculator=nm, n=nrow(dd), adjR2=round(summary(m)$adj.r.squared,3),
    resid_skew=round(skewness(r),2), resid_kurt=round(kurtosis(r)-3,2),
    BreuschPagan_p=signif(bp,2), SRI_nonlinearity_p=signif(p_nonlin,2), max_VIF=round(vmax,1))
  diaglist[[nm]]<-data.table(fit=fitted(m), res=rstandard(m), sri=dd[[v]], calc=nm)
}
tab<-rbindlist(rows); print(tab)
fwrite(tab, file.path(inter,"model_assumptions.csv"))

# The visual assumption checks (linearity, homogeneity, normality) are shown per
# calculator by performance::check_model() and car::qqPlot() in 98_check_model.R;
# this script only produces the quantitative summary table above.
cat("\nwritten model_assumptions.csv\nDONE\n")
