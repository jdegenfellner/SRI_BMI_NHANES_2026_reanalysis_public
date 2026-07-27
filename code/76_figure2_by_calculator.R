# =====================================================================
# 76_figure2_by_calculator.R  ->  figures/figure2_by_calculator.pdf
# Replicates Figure 2 of Degenfellner et al. (2026) -- predicted BMI across the
# SRI range by sex and race/ethnicity from the fully adjusted interaction model
# (the paper's model 4: log(BMI) ~ SRI x sex + SRI x ethnicity + covariates,
# survey-weighted) -- fitted once per SRI calculator, so the calculator-dependence
# of the sex/ethnic effect modification is visible side by side (ethnicity in
# rows, calculator in columns). Design-based 95% prediction bands.
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

covs <- c("RIDAGEYR","DMDEDUC2_factor","INDHHIN2_factor","OCD150_factor",
          "alcohol_consumption","smoking_status","LBXVIDMS","DR1TKCAL_winsorized",
          "Depression_score_PHQ_9","DMDMARTL_factor","activity_level")
eth_lv    <- c("Mexican American","Other Hispanic","Non-Hispanic White",
               "Non-Hispanic Black","Non-Hispanic Asian","Other Race - Including Multi-Racial")
eth_short <- c("Mexican Am.","Other Hisp.","NH White","NH Black","NH Asian","Other/Multi")

fitpred <- function(col){
  dd <- d[is.finite(get(col)) & is.finite(BMXBMI_winsorized) & !is.na(WTMEC4YR)]
  dd[, SRI := get(col)]
  des <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC4YR, nest=TRUE, data=dd)
  f <- as.formula(paste0("log(BMXBMI_winsorized) ~ SRI*RIDRETH3_factor + SRI*RIAGENDR_factor + ",
                         paste(covs, collapse=" + ")))
  m <- svyglm(f, design=des)
  # prediction grid: SRI 0-100 x sex x ethnicity; other covariates held at
  # reference (continuous at the sample mean, factors/characters at their mode)
  grid <- CJ(SRI=seq(0,100,by=2.5), .sex=levels(dd$RIAGENDR_factor), .eth=eth_lv, sorted=FALSE)
  grid[, RIAGENDR_factor := factor(.sex, levels=levels(dd$RIAGENDR_factor))]
  grid[, RIDRETH3_factor := factor(.eth, levels=levels(dd$RIDRETH3_factor))]
  for(v in covs){
    x <- dd[[v]]
    if(is.factor(x))         set(grid, j=v, value=factor(levels(x)[1], levels=levels(x)))
    else if(is.character(x)) set(grid, j=v, value=names(sort(table(x), decreasing=TRUE))[1])
    else                     set(grid, j=v, value=mean(x, na.rm=TRUE))
  }
  prd <- predict(m, newdata=as.data.frame(grid), type="link", se.fit=TRUE)
  pv  <- as.data.frame(prd); est <- pv[[1]]; se <- pv[[2]]
  grid[, `:=`(BMI=exp(est), lo=exp(est-1.96*se), hi=exp(est+1.96*se), calc=col)]
  setnames(grid, c("RIAGENDR_factor","RIDRETH3_factor"), c("Sex","Ethnicity"))
  grid[, .(SRI,Sex,Ethnicity,BMI,lo,hi,calc)]
}

res <- rbindlist(lapply(c("Published","GGIR","sleepreg"), fitpred))
res[, calc     := factor(calc, levels=c("Published","GGIR","sleepreg"))]
res[, Sex      := factor(Sex,  levels=c("Male","Female"))]
res[, Ethnicity:= factor(eth_short[match(as.character(Ethnicity), eth_lv)], levels=eth_short)]

female_col <- "#00BFC4"; male_col <- "#F8766D"
# fixed y-scale across all panels (as in the original Figure 2), tightly bounded so
# every row -- including the lower-BMI Non-Hispanic Asian row -- is visible while
# the band proportions match the published bootstrap figure
p <- ggplot(res, aes(SRI, BMI, color=Sex, fill=Sex)) +
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.18, colour=NA) +
  geom_line(linewidth=0.8) +
  facet_grid(Ethnicity ~ calc) +
  scale_color_manual(values=c(Male=male_col, Female=female_col)) +
  scale_fill_manual(values=c(Male=male_col, Female=female_col)) +
  coord_cartesian(ylim=c(19, 33)) +
  labs(x="Sleep Regularity Index", y="Predicted BMI", color="Sex", fill="Sex") +
  theme_bw(base_size=9) +
  theme(legend.position="top", panel.grid.minor=element_blank(),
        strip.background=element_rect(fill="grey92"), strip.text=element_text(size=8))
ggsave(file.path(figdir,"figure2_by_calculator.pdf"), p, width=8.2, height=9.2)
cat("wrote", file.path(figdir,"figure2_by_calculator.pdf"),"\n")
cat("BMI prediction range (within ylim clip):\n"); print(res[, .(minBMI=round(min(BMI),1), maxBMI=round(max(BMI),1)), by=calc])
