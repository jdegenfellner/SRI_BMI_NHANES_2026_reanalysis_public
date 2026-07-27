# =====================================================================
# performance::check_model() and car::qqPlot() for the fully adjusted model
# (Model 2 of Degenfellner et al. 2026, Table 3) re-fitted with each SRI
# calculator -- exactly the diagnostic displays used in the original paper.
# =====================================================================
suppressMessages({library(data.table); library(performance); library(see)
                  library(qqplotr); library(car); library(ggplot2)})
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())
source(file.path(base, if (dir.exists(file.path(base, "code"))) "code" else "01_recompute_SRI", "00_paths.R"))
# sdir provided by 00_paths.R

d  <- as.data.table(readRDS(file.path(origdata, "df_an1_imputed_MAIN.RDS")))
sri<- as.data.table(readRDS(file.path(inter,"sri_all_methods.rds")))[,.(SEQN,SRI_GGIR,SRI_sleepreg)]
d  <- merge(d, sri, by="SEQN", all.x=TRUE)
covs <- "RIDAGEYR + RIAGENDR_factor + RIDRETH3_factor + DMDEDUC2_factor + INDHHIN2_factor + OCD150_factor + alcohol_consumption + smoking_status + LBXVIDMS + DR1TKCAL_winsorized + Depression_score_PHQ_9 + DMDMARTL_factor + activity_level"
calcs <- c(Published="SRI_no_imputation_winsorized_rand", GGIR="SRI_GGIR", sleepreg="SRI_sleepreg")
# Build the formula as an object with its own environment (holding the data), so
# that check_model()/car refits resolve the model without needing `covs` in scope
# -- this keeps the script robust when sourced from run_all.R in a fresh env.
fitm <- function(v){
  dd <- d[is.finite(get(v))&is.finite(BMXBMI)]
  f  <- as.formula(paste0("log(BMXBMI) ~ ", v, " + ", covs)); environment(f) <- environment()
  lm(f, data=dd)
}

# ---- check_model() per calculator ----
for(nm in names(calcs)){
  m <- fitm(calcs[[nm]])
  p <- plot(check_model(m, verbose=FALSE))           # see:::plot.check_model -> patchwork
  ggsave(file.path(sdir, paste0("S_check_model_",nm,".pdf")), p, width=11, height=13, limitsize=FALSE)
  cat("check_model done:", nm, "\n")
}

# ---- car::qqPlot() 1x3 (base graphics) ----
pdf(file.path(sdir,"S_qqPlot_car.pdf"), width=12, height=4.5)
par(mfrow=c(1,3), mar=c(4,4,3,1))
for(nm in names(calcs)){
  dd <- d[is.finite(get(calcs[[nm]]))&is.finite(BMXBMI)]
  dd[, SRIvar := get(calcs[[nm]])]
  f <- as.formula(paste0("log(BMXBMI) ~ SRIvar + ", covs)); environment(f) <- environment()
  m <- lm(f, data=dd)
  car::qqPlot(m, envelope=list(level=0.95), id=FALSE, pch=20, col.lines="red",
              main=paste0("car::qqPlot - ", nm), ylab="studentized residuals", grid=FALSE)
}
dev.off()
cat("qqPlot (car) 2x2 written\nDONE\n")
