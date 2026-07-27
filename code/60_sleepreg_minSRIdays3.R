# =====================================================================
# sleepreg with minSRIdays = 3 (matching the original paper's >=3-day rule)
# instead of the package default 5. Re-uses the SWS files from part 2.
# Reports recovered n and the fully-adjusted BMI effect size (MF + CI).
# =====================================================================
suppressMessages({library(data.table); library(survey)})
options(width = 170)
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())
rdir  <- file.path(base,"01_recompute_SRI"); inter <- file.path(rdir,"intermediate")
source(file.path(rdir,"sleepreg_src","sleepreg_allFunctions.R"))
load(file.path(rdir,"sleepreg_src","sysdata.rda"))

swsdir <- file.path(inter,"sleepreg_binary_in","SWS_output")
out3   <- file.path(inter,"sleepreg_out_min3"); unlink(out3,recursive=TRUE); dir.create(out3)
invisible(capture.output(
  SRI_from_binary(binarydir=swsdir, alloutdir=out3, col.trans=1, col.timestamp=2,
                  wr.raster=FALSE, overwr=TRUE, minSRIdays=3)   # <-- the only change
))
sr3 <- fread(file.path(out3,"SRI.csv"))
sr3[, SEQN := as.numeric(gsub("\\D","",File))]
sr3 <- sr3[, .(SEQN, SRI_sleepreg3 = as.numeric(SRI), days3 = as.numeric(SRI_days))]
saveRDS(sr3, file.path(inter,"sri_sleepreg_min3.rds"))

# ---- modeling frame ----
d  <- as.data.table(readRDS(file.path(base,"00_original_paper_materials/data/df_an1_imputed_MAIN.RDS")))
an <- as.data.table(readRDS(file.path(base,"00_original_paper_materials/data/Combined_data_sets_2_waves_analytic_data_set.RDS")))
d  <- merge(d, an[,.(SEQN,WTMEC2YR,SDMVPSU,SDMVSTRA)], by="SEQN", all.x=TRUE)[, WTMEC4YR := WTMEC2YR/2]
d  <- merge(d, sr3[,.(SEQN,SRI_sleepreg3)], by="SEQN", all.x=TRUE)
# also bring sleepreg-5 for direct contrast
sr5 <- as.data.table(readRDS(file.path(inter,"sri_all_methods.rds")))[,.(SEQN,SRI_sleepreg)]
d  <- merge(d, sr5, by="SEQN", all.x=TRUE)

set.seed(3455); SRIsim <- replicate(100000, -100+200*sum(rbinom(1440,1,.5))/1440)
m_sri<-mean(SRIsim); s_sri<-sd(SRIsim)
winz <- function(x){ i<-which(x< -6.8); x[i]<-rnorm(length(i),m_sri,s_sri); x }
d[, SRI_sleepreg3_w := winz(SRI_sleepreg3)]
d[, SRI_sleepreg_w  := winz(SRI_sleepreg)]

covars <- paste("RIDAGEYR","RIAGENDR_factor","RIDRETH3_factor","DMDEDUC2_factor","INDHHIN2_factor",
                "OCD150_factor","alcohol_consumption","smoking_status","LBXVIDMS","DR1TKCAL_winsorized",
                "Depression_score_PHQ_9","DMDMARTL_factor","activity_level", sep=" + ")

fit <- function(svar){
  dd <- d[!is.na(get(svar)) & !is.na(BMXBMI) & !is.na(WTMEC4YR)]
  dd[, SRIz := scale(get(svar))[,1]]
  qb <- quantile(dd[[svar]],0:5/5,na.rm=TRUE); qb[1]<--Inf; qb[6]<-Inf
  dd[, SRIq := relevel(factor(cut(get(svar),qb,labels=paste0("Q",1:5),include.lowest=TRUE)),"Q1")]
  des <- svydesign(id=~SDMVPSU,strata=~SDMVSTRA,weights=~WTMEC4YR,nest=TRUE,data=dd)
  bc <- unname(summary(svyglm(as.formula(paste0("log(BMXBMI)~SRIz+",covars)),design=des))$coef["SRIz",])
  bq <- unname(summary(svyglm(as.formula(paste0("log(BMXBMI)~SRIq+",covars)),design=des))$coef["SRIqQ5",])
  data.table(n=nrow(dd),
    MF_perSD=round(exp(bc[1]),4), perSD_lo=round(exp(bc[1]-1.96*bc[2]),4), perSD_hi=round(exp(bc[1]+1.96*bc[2]),4),
    MF_Q5vQ1=round(exp(bq[1]),4), Q5_lo=round(exp(bq[1]-1.96*bq[2]),4), Q5_hi=round(exp(bq[1]+1.96*bq[2]),4))
}

cat("=== sleepreg recovered n by min-days rule (analytic 7085) ===\n")
cat("minSRIdays=5 (default): n =", sum(!is.na(d$SRI_sleepreg)),
    " | minSRIdays=3: n =", sum(!is.na(d$SRI_sleepreg3)), "\n\n")
cat("=== fully-adjusted BMI effect (MF + 95% CI) ===\n")
res <- rbind(
  cbind(method="sleepreg (min5)", fit("SRI_sleepreg_w")),
  cbind(method="sleepreg (min3)", fit("SRI_sleepreg3_w")))
print(res)
# reference rows from the main comparison
mc <- as.data.table(readRDS(file.path(inter,"model_comparison.rds")))
cat("\n=== reference (from main comparison) ===\n")
print(mc[method %in% c("orig","GGIR"), .(method, n, MF_perSD, perSD_lo, perSD_hi, MF_Q5vQ1, Q5_lo, Q5_hi)])
fwrite(res, file.path(inter,"sleepreg_min3_models.csv"))
cat("\nDONE\n")
