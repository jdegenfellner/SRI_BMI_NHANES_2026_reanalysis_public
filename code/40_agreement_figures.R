# =====================================================================
# PART 4: Agreement between SRI methods — figures + quintile reshuffling.
# Differences are judged on the DISTRIBUTION OF THE DIFFERENCE
# (Bland-Altman), never on overlap of separate CIs.
# =====================================================================
suppressMessages({library(data.table); library(ggplot2)})
options(width = 160)

base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())
source(file.path(base, if (dir.exists(file.path(base, "code"))) "code" else "01_recompute_SRI", "00_paths.R"))
# figdir provided by 00_paths.R

all <- as.data.table(readRDS(file.path(inter, "sri_all_methods.rds")))
# analytic frame only (the 7085)
a <- all[!is.na(SRI_orig)]
methcols <- c(Hand="SRI_hand", sleepreg="SRI_sleepreg", GGIR="SRI_GGIR", Original="SRI_orig")

# ---- attrition by method (how many get an SRI at all) ----
cat("=== valid SRI per method (analytic n=", nrow(a), ") ===\n", sep="")
for (nm in names(methcols)) cat(sprintf("  %-9s %d valid (%d NA)\n", nm,
     sum(!is.na(a[[methcols[[nm]]]])), sum(is.na(a[[methcols[[nm]]]]))))

# ---- pairwise correlation matrix ----
cc <- a[, ..methcols]; setnames(cc, names(methcols))
cat("\n=== Pearson correlations ===\n"); print(round(cor(cc, use="pairwise.complete.obs"), 3))
cat("\n=== Spearman correlations ===\n"); print(round(cor(cc, use="pairwise.complete.obs", method="spearman"), 3))

# ---- Bland-Altman of each package vs the hand calculation ----
ba_plot <- function(x, y, xl, yl){
  d <- data.table(m=(x+y)/2, diff=y-x)[is.finite(m) & is.finite(diff)]
  mu <- mean(d$diff); s <- sd(d$diff); lo <- mu-1.96*s; hi <- mu+1.96*s
  ggplot(d, aes(m, diff)) + geom_point(alpha=.15, size=.5) +
    geom_hline(yintercept=mu, color="blue") +
    geom_hline(yintercept=c(lo,hi), color="red", linetype="dashed") +
    labs(x=paste0("Mean of ", xl, " & ", yl, " SRI"),
         y=paste0(yl, " - ", xl),
         subtitle=sprintf("bias=%.2f, 95%% LoA [%.1f, %.1f]", mu, lo, hi)) +
    theme_minimal()
}
for (nm in c("sleepreg","GGIR")){
  p <- ba_plot(a$SRI_hand, a[[methcols[[nm]]]], "Hand", nm)
  ggsave(file.path(figdir, paste0("BlandAltman_", nm, "_vs_Hand.pdf")), p, width=6, height=4.5)
}
# scatter matrix-ish: each package vs hand
for (nm in c("sleepreg","GGIR","Original")){
  d <- a[is.finite(SRI_hand) & is.finite(get(methcols[[nm]]))]
  p <- ggplot(d, aes(SRI_hand, get(methcols[[nm]]))) + geom_point(alpha=.15,size=.5) +
       geom_abline(slope=1,intercept=0,color="red",linetype="dashed") +
       geom_smooth(method="lm",se=FALSE,color="blue") +
       labs(x="Hand SRI", y=paste0(nm," SRI")) + theme_minimal() + coord_equal()
  ggsave(file.path(figdir, paste0("scatter_", nm, "_vs_Hand.pdf")), p, width=4.5, height=4.5)
}

# ---- quintile reshuffling: hand vs each package ----
qcut <- function(x){ qb<-quantile(x,0:5/5,na.rm=TRUE); qb[1]<--Inf; qb[6]<-Inf
                     cut(x,qb,labels=paste0("Q",1:5),include.lowest=TRUE) }
a[, q_hand := qcut(SRI_hand)]
cat("\n=== quintile reshuffling vs Hand (share staying in same quintile) ===\n")
for (nm in c("sleepreg","GGIR")){
  a[, qx := qcut(get(methcols[[nm]]))]
  ok <- a[!is.na(q_hand)&!is.na(qx)]
  same <- mean(ok$q_hand==ok$qx)
  cat(sprintf("  %-9s: %.1f%% same quintile as Hand (n=%d); off-by>=2: %.1f%%\n",
      nm, 100*same, nrow(ok),
      100*mean(abs(as.integer(ok$q_hand)-as.integer(ok$qx))>=2)))
  print(table(Hand=ok$q_hand, pkg=ok$qx))
}
cat("\nfigures written to", figdir, "\nDONE part 4\n")
