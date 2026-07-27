# =====================================================================
# 75_subgroup_diff.R  ->  figures/forest_subgroup_diff.pdf
# Paired-bootstrap DIFFERENCE in Model-2 MF (package - published), within each
# recording-completeness subgroup (seven consecutive days vs a day gap). The same
# participants are resampled for both calculators in each replicate. Point and
# 95% interval; the dashed line at 0 marks no difference (a reference, not a
# threshold; no stars). Companion panel B to forest_subgroup.pdf.
# =====================================================================
suppressMessages({library(data.table); library(survey); library(ggplot2)})
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())
source(file.path(base, if (dir.exists(file.path(base, "code"))) "code" else "01_recompute_SRI", "00_paths.R"))

d  <- as.data.table(readRDS(insrc("df_an1_imputed_MAIN.RDS")))
an <- as.data.table(readRDS(insrc("Combined_data_sets_2_waves_analytic_data_set.RDS")))
d  <- merge(d, an[,.(SEQN,WTMEC2YR,SDMVPSU,SDMVSTRA)], by="SEQN", all.x=TRUE)[, WTMEC4YR := WTMEC2YR/2]
sri<- as.data.table(readRDS(insrc("sri_all_methods.rds")))
d  <- merge(d, sri[,.(SEQN,SRI_GGIR,SRI_sleepreg)], by="SEQN", all.x=TRUE)
nd <- as.data.table(read.csv(insrc("participant_ndays.csv")))
d  <- merge(d, nd[,.(SEQN,n_days)], by="SEQN", all.x=TRUE)
set.seed(3455); SS<-replicate(100000,-100+200*sum(rbinom(1440,1,.5))/1440); m_s<-mean(SS); s_s<-sd(SS)
winz<-function(x){i<-which(x< -6.8);x[i]<-rnorm(length(i),m_s,s_s);x}
d[, `:=`(Published=SRI_no_imputation_winsorized_rand, GGIR=winz(SRI_GGIR), sleepreg=winz(SRI_sleepreg))]
cov_full <- c("RIDAGEYR","RIAGENDR_factor","RIDRETH3_factor","DMDEDUC2_factor","INDHHIN2_factor",
  "OCD150_factor","alcohol_consumption","smoking_status","LBXVIDMS","DR1TKCAL_winsorized",
  "Depression_score_PHQ_9","DMDMARTL_factor","activity_level")
rhs<-paste(cov_full,collapse=" + ")
delta<-function(x){q<-quantile(x,0:5/5,na.rm=TRUE);q[1]<-q[1]-1e-9
  m<-tapply(x,cut(x,q,include.lowest=TRUE,labels=1:5),median,na.rm=TRUE);as.numeric(m[5]-m[1])}
dsg<-function(dd) svydesign(id=~SDMVPSU,strata=~SDMVSTRA,weights=~WTMEC4YR,nest=TRUE,data=dd)
mf2<-function(dd,col,D){b<-coef(svyglm(as.formula(paste0("log(BMXBMI) ~ ",col," + ",rhs)),design=dsg(dd)))[col];exp(unname(b)*D)}
pair<-function(sub,colA,nb=500){
  dd<-d[sub & is.finite(get(colA))&is.finite(Published)&!is.na(BMXBMI)&!is.na(WTMEC4YR)]
  DA<-delta(dd[[colA]]);DB<-delta(dd$Published);n<-nrow(dd)
  pt<-mf2(dd,colA,DA)-mf2(dd,"Published",DB);v<-numeric(nb)
  for(i in 1:nb){bs<-dd[sample(n,replace=TRUE)];v[i]<-tryCatch(mf2(bs,colA,DA)-mf2(bs,"Published",DB),error=function(e)NA)}
  ci<-quantile(v,c(.025,.975),na.rm=TRUE); data.table(calc=colA,diff=pt,lo=ci[1],hi=ci[2]) }
set.seed(11)
grps <- list("Seven consecutive days (60%)"=quote(d$n_days==7), "Day gap (40%)"=quote(d$n_days<7))
res <- rbindlist(lapply(names(grps), function(gn){ sub<-eval(grps[[gn]])
  rbindlist(lapply(c("GGIR","sleepreg"), function(cc){ r<-pair(sub,cc); r[,subgroup:=gn]; r })) }))
res[, subgroup:=factor(subgroup,levels=names(grps))]
res[, calc:=factor(calc,levels=c("sleepreg","GGIR"))]

p<-ggplot(res, aes(diff, calc))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey55")+
  geom_errorbar(aes(xmin=lo,xmax=hi), orientation="y", width=.18, colour="grey30")+
  geom_point(aes(colour=subgroup), size=2.6)+
  geom_text(aes(label=sprintf("%+.3f",diff)), vjust=-1.0, size=3)+
  facet_wrap(~subgroup, ncol=1)+
  scale_colour_manual(values=c("#2c7fb8","#d95f0e"), guide="none")+
  labs(x="Difference in BMI multiplication factor vs the published SRI\n(paired bootstrap, 95% interval)", y=NULL)+
  theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),
        strip.background=element_rect(fill="grey92"), strip.text=element_text(face="bold"))
ggsave(file.path(figdir,"forest_subgroup_diff.pdf"), p, width=6.6, height=4.2)
cat("wrote", file.path(figdir,"forest_subgroup_diff.pdf"),"\n"); print(res[,.(subgroup,calc,diff=round(diff,3),lo=round(lo,3),hi=round(hi,3))])
