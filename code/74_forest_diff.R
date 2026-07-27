# =====================================================================
# 74_forest_diff.R  ->  figures/forest_diff.pdf
# Forest of the DIFFERENCE in MF vs the published SRI (paired bootstrap),
# for Model 2 and Model 3 (women). Point + 95% interval; dashed line at 0 is a
# reference for "no difference", NOT a significance threshold (no stars).
# =====================================================================
suppressMessages({library(data.table); library(survey); library(ggplot2)})
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())
source(file.path(base, if (dir.exists(file.path(base, "code"))) "code" else "01_recompute_SRI", "00_paths.R"))

d  <- as.data.table(readRDS(insrc("df_an1_imputed_MAIN.RDS")))
an <- as.data.table(readRDS(insrc("Combined_data_sets_2_waves_analytic_data_set.RDS")))
d  <- merge(d, an[,.(SEQN,WTMEC2YR,SDMVPSU,SDMVSTRA)], by="SEQN", all.x=TRUE)[, WTMEC4YR := WTMEC2YR/2]
sri<- as.data.table(readRDS(insrc("sri_all_methods.rds")))
d  <- merge(d, sri[,.(SEQN,SRI_GGIR,SRI_sleepreg)], by="SEQN", all.x=TRUE)
set.seed(3455); SS<-replicate(100000,-100+200*sum(rbinom(1440,1,.5))/1440); m_s<-mean(SS); s_s<-sd(SS)
winz<-function(x){i<-which(x< -6.8);x[i]<-rnorm(length(i),m_s,s_s);x}
d[, `:=`(Published=SRI_no_imputation_winsorized_rand, GGIR=winz(SRI_GGIR), sleepreg=winz(SRI_sleepreg))]
cov_full <- c("RIDAGEYR","RIAGENDR_factor","RIDRETH3_factor","DMDEDUC2_factor","INDHHIN2_factor",
  "OCD150_factor","alcohol_consumption","smoking_status","LBXVIDMS","DR1TKCAL_winsorized",
  "Depression_score_PHQ_9","DMDMARTL_factor","activity_level")
rhs<-paste(cov_full,collapse=" + "); rhs_ns<-paste(setdiff(cov_full,"RIAGENDR_factor"),collapse=" + ")
delta<-function(x){q<-quantile(x,0:5/5,na.rm=TRUE);q[1]<-q[1]-1e-9
  m<-tapply(x,cut(x,q,include.lowest=TRUE,labels=1:5),median,na.rm=TRUE);as.numeric(m[5]-m[1])}
dsg<-function(dd) svydesign(id=~SDMVPSU,strata=~SDMVSTRA,weights=~WTMEC4YR,nest=TRUE,data=dd)
mf2<-function(dd,col,D){b<-coef(svyglm(as.formula(paste0("log(BMXBMI) ~ ",col," + ",rhs)),design=dsg(dd)))[col];exp(unname(b)*D)}
mfw<-function(dd,col,D){cf<-coef(svyglm(as.formula(paste0("log(BMXBMI) ~ RIAGENDR_factor*",col," + ",rhs_ns)),design=dsg(dd)))
  ix<-grep(paste0(":",col,"$|^",col,":"),names(cf),value=TRUE);bfem<-cf[col]+(if(length(ix))cf[ix[1]] else 0);exp(unname(bfem)*D)}
pair<-function(colA,colB,fn,nb=500){
  dd<-d[is.finite(get(colA))&is.finite(get(colB))&!is.na(BMXBMI)&!is.na(WTMEC4YR)]
  DA<-delta(dd[[colA]]);DB<-delta(dd[[colB]]);n<-nrow(dd);pt<-fn(dd,colA,DA)-fn(dd,colB,DB);v<-numeric(nb)
  for(i in 1:nb){bs<-dd[sample(n,replace=TRUE)];v[i]<-tryCatch(fn(bs,colA,DA)-fn(bs,colB,DB),error=function(e)NA)}
  ci<-quantile(v,c(.025,.975),na.rm=TRUE); data.table(calc=colA, diff=pt, lo=ci[1], hi=ci[2]) }
set.seed(11)
res<-rbindlist(list(
  pair("GGIR","Published",mf2)[,panel:="Fully adjusted (Model 2)"],
  pair("sleepreg","Published",mf2)[,panel:="Fully adjusted (Model 2)"],
  pair("GGIR","Published",mfw)[,panel:="Women (Model 3)"],
  pair("sleepreg","Published",mfw)[,panel:="Women (Model 3)"]))
res[, calc:=factor(calc, levels=c("sleepreg","GGIR"))]
res[, panel:=factor(panel, levels=c("Fully adjusted (Model 2)","Women (Model 3)"))]

p<-ggplot(res, aes(diff, calc))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey55")+
  geom_errorbar(aes(xmin=lo,xmax=hi), orientation="y", width=.18, colour="grey30")+
  geom_point(size=2.6, colour="#d95f0e")+
  geom_text(aes(label=sprintf("%+.3f",diff)), vjust=-1.0, size=3)+
  facet_wrap(~panel)+
  labs(x="Difference in BMI multiplication factor vs the published SRI (paired bootstrap, 95% interval)", y=NULL)+
  theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),
        strip.background=element_rect(fill="grey92"), strip.text=element_text(face="bold"))
ggsave(file.path(figdir,"forest_diff.pdf"), p, width=8.6, height=2.4)
cat("wrote", file.path(figdir,"forest_diff.pdf"),"\n"); print(res[,.(panel,calc,diff=round(diff,3),lo=round(lo,3),hi=round(hi,3))])
