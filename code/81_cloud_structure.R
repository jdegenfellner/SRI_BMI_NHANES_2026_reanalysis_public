# =====================================================================
# Why does the Hand-vs-package SRI scatter show several clouds / slopes?
# Colour the scatter by candidate per-participant features to find the
# latent variable that separates the clouds.
# =====================================================================
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())
source(file.path(base, if (dir.exists(file.path(base, "code"))) "code" else "01_recompute_SRI", "00_paths.R"))

all <- as.data.table(readRDS(insrc("sri_all_methods.rds")))[!is.na(SRI_orig)]
sr  <- as.data.table(readRDS(insrc("sri_sleepreg.rds")))           # sleepreg_days
# per-participant features (valid days, missingness, day-pairing structure);
# precomputed by 11_participant_features.R so the minute series is not needed
feat <- as.data.table(read.csv(insrc("participant_ndays.csv")))[, .(SEQN, n_days, miss, n_adj, n_gaps)]

d <- Reduce(function(a,b) merge(a,b,by="SEQN",all.x=TRUE),
            list(all[,.(SEQN,SDDSRVYR,SRI_orig,SRI_GGIR,SRI_sleepreg)], feat, sr[,.(SEQN,sleepreg_days)]))
d[, wave := factor(SDDSRVYR, labels=c("2011-12","2013-14"))]
d[, has_gap := factor(n_gaps>0, labels=c("all days consecutive","has a day gap"))]
cat("=== feature summary ===\n")
print(d[, .(n=.N, med_days=median(n_days), med_miss=round(median(miss),3),
            share_with_gap=round(mean(n_gaps>0),3))])
cat("\n=== n_days distribution ===\n"); print(table(d$n_days))
cat("=== n_gaps (non-consecutive day boundaries) ===\n"); print(table(d$n_gaps))

# ---- residual analysis: SRI_sleepreg ~ SRI_orig, what predicts the residual?
d2 <- d[is.finite(SRI_sleepreg)]
d2[, resid := SRI_sleepreg - SRI_orig]
cat("\n=== correlation of (sleepreg - hand) gap with candidate features ===\n")
for(v in c("n_days","miss","n_gaps","sleepreg_days")){
  cat(sprintf("  %-14s : r = %+.3f\n", v, cor(d2$resid, d2[[v]], use="complete.obs")))
}
cat("mean gap by has_gap:\n"); print(d2[, .(mean_gap=round(mean(resid),2), n=.N), by=has_gap])
cat("mean gap by n_days:\n"); print(d2[, .(mean_gap=round(mean(resid),2), n=.N), by=n_days][order(n_days)])

# ---- coloured scatters (Hand vs sleepreg) by each candidate ----
base_pts <- function(colvar, lab, cont=TRUE){
  p <- ggplot(d[is.finite(SRI_sleepreg)], aes(SRI_orig, SRI_sleepreg, colour=.data[[colvar]])) +
    geom_point(alpha=.25, size=.5) +
    geom_abline(slope=1,intercept=0,linetype="dashed",colour="black") +
    labs(x="Published SRI", y="sleepreg SRI", colour=lab) + theme_minimal(base_size=10)
  if(cont) p <- p + scale_colour_viridis_c() else p <- p + scale_colour_brewer(palette="Set1")
  p
}
pg <- (base_pts("n_days","valid days",TRUE) | base_pts("miss","missingness",TRUE)) /
      (base_pts("has_gap","day gap?",FALSE) | base_pts("wave","NHANES wave",FALSE))
ggsave(file.path(figdir,"cloud_structure_Hand_vs_sleepreg.pdf"), pg, width=11, height=9)
# same for GGIR
base_pts2 <- function(colvar, lab, cont=TRUE){
  p <- ggplot(d, aes(SRI_orig, SRI_GGIR, colour=.data[[colvar]])) +
    geom_point(alpha=.25,size=.5) + geom_abline(slope=1,intercept=0,linetype="dashed",colour="black") +
    labs(x="Published SRI", y="GGIR SRI", colour=lab) + theme_minimal(base_size=10)
  if(cont) p+scale_colour_viridis_c() else p+scale_colour_brewer(palette="Set1")
}
pg2 <- (base_pts2("n_days","valid days",TRUE) | base_pts2("miss","missingness",TRUE)) /
       (base_pts2("has_gap","day gap?",FALSE) | base_pts2("wave","NHANES wave",FALSE))
ggsave(file.path(figdir,"cloud_structure_Hand_vs_GGIR.pdf"), pg2, width=11, height=9)
cat("\nfigures written\nDONE\n")
