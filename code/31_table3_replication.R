# =====================================================================
# DIRECT replication of the published Table 3 method, for every SRI
# calculator, so each cell is comparable to Degenfellner et al. (2026).
#
# Table-3 recipe (from original 4_Tables.R):
#  * Models (survey-weighted, CONTINUOUS winsorized SRI):
#     M1 = age + sex ; M2 = full (13 covars, no interaction) ;
#     M3 = full + sex x SRI ; M4 = full + sex x SRI + ethnicity x SRI
#  * Quintiles via ntile(SRI,5); take MEDIAN SRI per quintile.
#  * MF Q1->Qk = exp( beta_SRI * (median_Qk - median_Q1) )   [linear in SRI;
#    covariates cancel for non-interaction models; for interaction models
#    beta_SRI is the sex/ethnicity-specific slope].
#  * 95% CI by bootstrap (500 resamples, refit, medians fixed) -- as original.
# Validation target (orig SRI): M1 Q1->Q5 = 0.9225 ; M2 = 0.9503 ;
#    M3 female = 0.92 , male = 0.981.
# =====================================================================
suppressMessages({library(data.table); library(survey)})
options(width=170)
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())
source(file.path(base, if (dir.exists(file.path(base, "code"))) "code" else "01_recompute_SRI", "00_paths.R"))

d  <- as.data.table(readRDS(file.path(origdata, "df_an1_imputed_MAIN.RDS")))
an <- as.data.table(readRDS(file.path(origdata, "Combined_data_sets_2_waves_analytic_data_set.RDS")))
d  <- merge(d, an[,.(SEQN,WTMEC2YR,SDMVPSU,SDMVSTRA)], by="SEQN", all.x=TRUE)[, WTMEC4YR := WTMEC2YR/2]
sri <- as.data.table(readRDS(file.path(inter,"sri_all_methods.rds")))
d  <- merge(d, sri[,.(SEQN,SRI_hand,SRI_GGIR,SRI_sleepreg)], by="SEQN", all.x=TRUE)
sr3 <- as.data.table(readRDS(file.path(inter,"sri_sleepreg_min3.rds")))
d  <- merge(d, sr3[,.(SEQN,SRI_sleepreg3)], by="SEQN", all.x=TRUE)

set.seed(3455); SRIsim <- replicate(100000,-100+200*sum(rbinom(1440,1,.5))/1440)
m_s<-mean(SRIsim); s_s<-sd(SRIsim)
winz <- function(x){ i<-which(x< -6.8); x[i]<-rnorm(length(i),m_s,s_s); x }
d[, `:=`(orig=SRI_no_imputation_winsorized_rand, hand=winz(SRI_hand),
         GGIR=winz(SRI_GGIR), sleepreg5=winz(SRI_sleepreg), sleepreg3=winz(SRI_sleepreg3))]

cov_full <- c("RIDAGEYR","RIAGENDR_factor","RIDRETH3_factor","DMDEDUC2_factor","INDHHIN2_factor",
              "OCD150_factor","alcohol_consumption","smoking_status","LBXVIDMS","DR1TKCAL_winsorized",
              "Depression_score_PHQ_9","DMDMARTL_factor","activity_level")
cov_full_nosex <- setdiff(cov_full,"RIAGENDR_factor")

beta_SRI <- function(des, rhs) {
  m <- svyglm(as.formula(paste0("log(BMXBMI) ~ SRI + ", rhs)), design=des)
  coef(m)
}

# one calculator -> point MFs for M1,M2,M3(sex)
mf_for <- function(dd, svarvals){
  dd <- copy(dd); dd[, SRI := svarvals]
  dd <- dd[is.finite(SRI) & !is.na(BMXBMI) & !is.na(WTMEC4YR)]
  med <- dd[, median(SRI), by=.(q=cut(SRI,quantile(SRI,0:5/5,na.rm=TRUE),include.lowest=TRUE,labels=1:5))][order(q),V1]
  dSRI <- med[5]-med[1]
  des <- svydesign(id=~SDMVPSU,strata=~SDMVSTRA,weights=~WTMEC4YR,nest=TRUE,data=dd)
  b1 <- beta_SRI(des, "RIDAGEYR + RIAGENDR_factor")["SRI"]
  b2 <- beta_SRI(des, paste(cov_full,collapse=" + "))["SRI"]
  c3 <- beta_SRI(des, paste0("RIAGENDR_factor*SRI - SRI + SRI + ", paste(cov_full_nosex,collapse=" + ")))
  # sex-specific slopes from interaction model
  ix <- grep("RIAGENDR_factor.*:SRI|SRI:RIAGENDR_factor", names(c3), value=TRUE)
  bref <- c3["SRI"]; bother <- bref + (if(length(ix)) c3[ix[1]] else 0)
  list(n=nrow(dd), med=med, dSRI=dSRI,
       M1=exp(b1*dSRI), M2=exp(b2*dSRI),
       M3_ref=exp(bref*dSRI), M3_other=exp(bother*dSRI),
       reflev=levels(dd$RIAGENDR_factor)[1], otherlev=levels(dd$RIAGENDR_factor)[2])
}

# bootstrap CI (refit; medians fixed from full data) -- like original (n_boot=500)
boot_mf <- function(dd, svarvals, nb=500){
  dd <- copy(dd); dd[, SRI := svarvals]
  dd <- dd[is.finite(SRI) & !is.na(BMXBMI) & !is.na(WTMEC4YR)]
  med <- dd[, median(SRI), by=.(q=cut(SRI,quantile(SRI,0:5/5,na.rm=TRUE),include.lowest=TRUE,labels=1:5))][order(q),V1]
  dSRI <- med[5]-med[1]; n <- nrow(dd)
  M <- matrix(NA,nb,4, dimnames=list(NULL,c("M1","M2","M3_ref","M3_other")))
  for(i in 1:nb){
    bs <- dd[sample(n,replace=TRUE)]
    des <- svydesign(id=~SDMVPSU,strata=~SDMVSTRA,weights=~WTMEC4YR,nest=TRUE,data=bs)
    b1 <- tryCatch(beta_SRI(des,"RIDAGEYR + RIAGENDR_factor")["SRI"],error=function(e)NA)
    b2 <- tryCatch(beta_SRI(des,paste(cov_full,collapse=" + "))["SRI"],error=function(e)NA)
    c3 <- tryCatch(beta_SRI(des,paste0("RIAGENDR_factor*SRI - SRI + SRI + ",paste(cov_full_nosex,collapse=" + "))),error=function(e)NULL)
    if(!is.null(c3)){ ix<-grep("RIAGENDR_factor.*:SRI|SRI:RIAGENDR_factor",names(c3),value=TRUE)
      br<-c3["SRI"]; bo<-br+(if(length(ix))c3[ix[1]] else 0)
      M[i,]<-c(exp(b1*dSRI),exp(b2*dSRI),exp(br*dSRI),exp(bo*dSRI))
    } else M[i,]<-c(exp(b1*dSRI),exp(b2*dSRI),NA,NA)
  }
  apply(M,2,function(x) quantile(x,c(.025,.975),na.rm=TRUE))
}

calcs <- c("orig","hand","GGIR","sleepreg5","sleepreg3")
fmt <- function(p,lo,hi) sprintf("%.3f (%.3f; %.3f)",p,lo,hi)
set.seed(11)
rows <- list()
for(cc in calcs){
  pt <- mf_for(d, d[[cc]])
  ci <- boot_mf(d, d[[cc]], nb=500)
  rows[[cc]] <- data.table(
    calculator=cc, n=pt$n,
    M1_age_sex          = fmt(pt$M1,    ci["2.5%","M1"],     ci["97.5%","M1"]),
    M2_fully_adj        = fmt(pt$M2,    ci["2.5%","M2"],     ci["97.5%","M2"]),
    M3_sex_ref          = paste0(pt$reflev,": ",  fmt(pt$M3_ref,  ci["2.5%","M3_ref"],  ci["97.5%","M3_ref"])),
    M3_sex_other        = paste0(pt$otherlev,": ",fmt(pt$M3_other,ci["2.5%","M3_other"],ci["97.5%","M3_other"])))
  cat("done:",cc,"\n")
}
tab <- rbindlist(rows)
fwrite(tab, file.path(inter,"table3_replication.csv"))
saveRDS(tab, file.path(inter,"table3_replication.rds"))
cat("\n================= Table-3-style MF Q1->Q5 by SRI calculator =================\n")
for(i in 1:nrow(tab)){
  cat(sprintf("\n%-11s (n=%d)\n", tab$calculator[i], tab$n[i]))
  cat("  M1 age+sex  :", tab$M1_age_sex[i], "\n")
  cat("  M2 full adj :", tab$M2_fully_adj[i], "\n")
  cat("  M3 sex-split:", tab$M3_sex_ref[i], "|", tab$M3_sex_other[i], "\n")
}
cat("\nValidation target (published, orig): M1=0.9225  M2=0.9503  M3 Female=0.92 Male=0.981\n")
cat("DONE\n")
