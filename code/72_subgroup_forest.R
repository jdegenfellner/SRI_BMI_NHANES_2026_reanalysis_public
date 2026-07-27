# =====================================================================
# 72_subgroup_forest.R  ->  figures/forest_subgroup.pdf  (manuscript Fig.)
# Model-2 MF (Q1->Q5) per calculator, split by recording completeness.
# Shows the SRI-BMI effect is calculator-invariant among 7-consecutive-day
# participants and diverges only in the gapped subgroup (measurement error).
# Runs from the shipped derived data.
# =====================================================================
suppressMessages({library(data.table); library(survey); library(ggplot2)})
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())
source(file.path(base, if (dir.exists(file.path(base, "code"))) "code" else "01_recompute_SRI", "00_paths.R"))

d  <- as.data.table(readRDS(insrc("df_an1_imputed_MAIN.RDS")))
an <- as.data.table(readRDS(insrc("Combined_data_sets_2_waves_analytic_data_set.RDS")))
d  <- merge(d, an[,.(SEQN,WTMEC2YR,SDMVPSU,SDMVSTRA)], by="SEQN", all.x=TRUE)[, WTMEC4YR := WTMEC2YR/2]
sri <- as.data.table(readRDS(insrc("sri_all_methods.rds")))
d  <- merge(d, sri[,.(SEQN,SRI_hand,SRI_GGIR,SRI_sleepreg)], by="SEQN", all.x=TRUE)
nd <- as.data.table(read.csv(insrc("participant_ndays.csv")))
d  <- merge(d, nd[,.(SEQN,n_days)], by="SEQN", all.x=TRUE)

set.seed(3455); SRIsim <- replicate(100000,-100+200*sum(rbinom(1440,1,.5))/1440)
m_s<-mean(SRIsim); s_s<-sd(SRIsim)
winz <- function(x){ i<-which(x< -6.8); x[i]<-rnorm(length(i),m_s,s_s); x }
d[, `:=`(orig=SRI_no_imputation_winsorized_rand,
         GGIR=winz(SRI_GGIR), sleepreg=winz(SRI_sleepreg))]
cov_full <- c("RIDAGEYR","RIAGENDR_factor","RIDRETH3_factor","DMDEDUC2_factor","INDHHIN2_factor",
              "OCD150_factor","alcohol_consumption","smoking_status","LBXVIDMS","DR1TKCAL_winsorized",
              "Depression_score_PHQ_9","DMDMARTL_factor","activity_level")

mf_sub <- function(svarvals, keep){
  dd <- copy(d); dd[, SRI := svarvals]; dd <- dd[keep & is.finite(SRI) & !is.na(BMXBMI) & !is.na(WTMEC4YR)]
  qn <- quantile(dd$SRI,0:5/5,na.rm=TRUE); qn[1]<-qn[1]-1e-9
  med <- dd[, median(SRI), by=.(q=cut(SRI,qn,include.lowest=TRUE,labels=1:5))][order(q),V1]
  dSRI <- med[5]-med[1]
  des <- svydesign(id=~SDMVPSU,strata=~SDMVSTRA,weights=~WTMEC4YR,nest=TRUE,data=dd)
  m <- svyglm(as.formula(paste0("log(BMXBMI) ~ SRI + ", paste(cov_full,collapse=" + "))), design=des)
  b<-unname(coef(m)["SRI"]); se<-unname(sqrt(vcov(m)["SRI","SRI"])); tc<-qt(.975,degf(des))
  data.table(MF=exp(b*dSRI), lo=exp((b-tc*se)*dSRI), hi=exp((b+tc*se)*dSRI), n=nrow(dd))
}
labs  <- c(orig="Published", GGIR="GGIR", sleepreg="sleepreg")
grps  <- list("Seven consecutive days (60%)" = quote(d$n_days==7),
              "Day gap (40%)"                = quote(d$n_days<7))
res <- rbindlist(lapply(names(grps), function(gn){
  keep <- eval(grps[[gn]])
  rbindlist(lapply(names(labs), function(cc){ r<-mf_sub(d[[cc]],keep); r[,`:=`(subgroup=gn,calc=labs[[cc]])]; r }))
}))
res[, subgroup := factor(subgroup, levels=names(grps))]
res[, calc := factor(calc, levels=rev(unname(labs)))]

p <- ggplot(res, aes(MF, calc)) +
  geom_vline(xintercept=1, linetype="dashed", colour="grey55") +
  geom_errorbar(aes(xmin=lo, xmax=hi), orientation="y", width=.22, colour="grey30") +
  geom_point(aes(colour=subgroup), size=2.6) +
  geom_text(aes(label=sprintf("%.3f", MF)), vjust=-1.0, size=3) +
  facet_wrap(~subgroup, ncol=1) +
  scale_colour_manual(values=c("#2c7fb8","#d95f0e"), guide="none") +
  labs(x="BMI multiplication factor, top vs bottom SRI quintile (Model 2, 95% CI)", y=NULL) +
  coord_cartesian(xlim=c(0.92, 1.04)) +
  theme_bw(base_size=11) + theme(panel.grid.minor=element_blank(),
        strip.background=element_rect(fill="grey92"), strip.text=element_text(face="bold"))
ggsave(file.path(figdir,"forest_subgroup.pdf"), p, width=6.6, height=4.2)
cat("wrote", file.path(figdir,"forest_subgroup.pdf"), "\n"); print(res[,.(subgroup,calc,MF=round(MF,3),n)])
