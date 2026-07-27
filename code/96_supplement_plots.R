# =====================================================================
# Supplementary figures for the commentary (mirroring the original paper's
# supplement, but ACROSS the three SRI calculators):
#  S-dist   : SRI distribution by calculator (overall + by NHANES cycle)
#  S-BA     : Bland-Altman of each package vs hand
#  S-diag   : model diagnostics (residuals) for the fully adjusted BMI model
#  S-pred   : predicted BMI across SRI by calculator
# (ggpairs 4-way, faceted-by-days and cloud-structure already exist.)
# =====================================================================
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())
source(file.path(base, if (dir.exists(file.path(base, "code"))) "code" else "01_recompute_SRI", "00_paths.R"))
# sdir provided by 00_paths.R

all <- as.data.table(readRDS(file.path(inter,"sri_all_methods.rds")))[!is.na(SRI_orig)]
long <- melt(all[,.(SEQN, SDDSRVYR, Published=SRI_orig, GGIR=SRI_GGIR, sleepreg=SRI_sleepreg)],
             id.vars=c("SEQN","SDDSRVYR"), variable.name="calculator", value.name="SRI")[is.finite(SRI)]
long[, cycle := factor(SDDSRVYR, labels=c("2011-12","2013-14"))]

# S-dist: density by calculator, overall + by cycle
p1 <- ggplot(long, aes(SRI, colour=calculator, fill=calculator)) +
  geom_density(alpha=.15, linewidth=.7) +
  labs(x="Sleep Regularity Index", y="density", title="A  SRI distribution by calculator (analytic sample)") +
  theme_minimal(base_size=11)+theme(plot.title=element_text(size=10))
p2 <- ggplot(long, aes(SRI, colour=calculator)) + geom_density(linewidth=.7) + facet_wrap(~cycle) +
  labs(x="Sleep Regularity Index", y="density", title="B  by NHANES cycle") +
  theme_minimal(base_size=11)+theme(plot.title=element_text(size=10))
ggsave(file.path(sdir,"S_SRI_distribution.pdf"), p1/p2, width=8.5, height=7)

# S-BA: Bland-Altman each package vs the published SRI
ba <- function(x,y,lab){ d<-data.table(m=(x+y)/2,diff=y-x)[is.finite(m)&is.finite(diff)]
  mu<-mean(d$diff);s<-sd(d$diff)
  ggplot(d,aes(m,diff))+geom_point(alpha=.12,size=.4)+
   geom_hline(yintercept=mu,colour="blue")+geom_hline(yintercept=c(mu-1.96*s,mu+1.96*s),colour="red",linetype="dashed")+
   labs(x="mean of the two SRIs",y=paste0(lab," - Published"),
        subtitle=sprintf("bias %.1f, 95%% LoA [%.1f, %.1f]",mu,mu-1.96*s,mu+1.96*s))+theme_minimal(base_size=10) }
ggsave(file.path(sdir,"S_BlandAltman.pdf"),
       ba(all$SRI_orig,all$SRI_GGIR,"GGIR") | ba(all$SRI_orig,all$SRI_sleepreg,"sleepreg"),
       width=9,height=4.3)

# S-diag + S-pred: fully adjusted BMI model per calculator (unweighted lm, for residual diagnostics)
d  <- as.data.table(readRDS(file.path(origdata, "df_an1_imputed_MAIN.RDS")))
sri<- as.data.table(readRDS(file.path(inter,"sri_all_methods.rds")))[,.(SEQN,SRI_GGIR,SRI_sleepreg)]
d  <- merge(d, sri, by="SEQN", all.x=TRUE)
covs <- "RIDAGEYR + RIAGENDR_factor + RIDRETH3_factor + DMDEDUC2_factor + INDHHIN2_factor + OCD150_factor + alcohol_consumption + smoking_status + LBXVIDMS + DR1TKCAL_winsorized + Depression_score_PHQ_9 + DMDMARTL_factor + activity_level"
diag_list <- list(); pred_list <- list()
for(nm in c("Published=SRI_no_imputation_winsorized_rand","GGIR=SRI_GGIR","sleepreg=SRI_sleepreg")){
  lab<-sub("=.*","",nm); v<-sub(".*=","",nm)
  dd<-d[is.finite(get(v))&is.finite(BMXBMI)]
  m<-lm(as.formula(paste0("log(BMXBMI) ~ ",v," + ",covs)), data=dd)
  rd<-data.table(fitted=fitted(m), resid=rstandard(m), calc=lab)
  diag_list[[lab]]<-rd
  gx<-seq(quantile(dd[[v]],.02,na.rm=TRUE),quantile(dd[[v]],.98,na.rm=TRUE),length=50)
  nd<-dd[rep(1,50)]; nd[[v]]<-gx
  for(cc in c("RIDAGEYR")) nd[[cc]]<-mean(dd[[cc]],na.rm=TRUE)
  pr<-exp(predict(m,newdata=nd)); pred_list[[lab]]<-data.table(SRI=gx, BMI=pr, calc=lab)
}
pred<-rbindlist(pred_list)
# (QQ + residual diagnostics now come from car::qqPlot and performance::check_model in 98_check_model.R)
ppred <- ggplot(pred, aes(SRI, BMI, colour=calc)) + geom_line(linewidth=1) +
  labs(title="Predicted BMI across SRI by calculator (fully adjusted, covariates at reference)",
       x="Sleep Regularity Index", y="predicted BMI", colour="calculator")+theme_minimal(base_size=10)+theme(plot.title=element_text(size=10))
ggsave(file.path(sdir,"S_predicted_BMI.pdf"), ppred, width=8, height=4.5)
cat("supplement figures written to", sdir, "\n")
cat("SRI medians by calculator x cycle:\n")
print(long[,.(median=round(median(SRI),1),n=.N),by=.(calculator,cycle)][order(calculator,cycle)])
cat("DONE\n")
