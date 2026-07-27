# =====================================================================
# Commentary figures: (1) forest plot of effect sizes (MF + 95% CI) across
# the SRI calculators; (2) agreement scatter of package SRIs vs hand.
# =====================================================================
suppressMessages({library(data.table); library(ggplot2)})
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())
source(file.path(base, if (dir.exists(file.path(base, "code"))) "code" else "01_recompute_SRI", "00_paths.R"))
# figdir provided by 00_paths.R

mc  <- as.data.table(readRDS(file.path(inter,"model_comparison.rds")))
m3  <- fread(file.path(inter,"sleepreg_min3_models.csv"))   # has method, MF_perSD, perSD_lo/hi, MF_Q5vQ1, Q5_lo/hi

lab <- c(orig="Original (hand, published)", hand="Hand (reproduction)",
         GGIR="GGIR", sleepreg="sleepreg (min 5 d)")
dd <- rbindlist(list(
  mc[, .(method=lab[method], perSD=MF_perSD, pl=perSD_lo, ph=perSD_hi, q=MF_Q5vQ1, ql=Q5_lo, qh=Q5_hi)],
  data.table(method="sleepreg (min 3 d)",
             perSD=m3[method=="sleepreg (min3)",MF_perSD], pl=m3[method=="sleepreg (min3)",perSD_lo], ph=m3[method=="sleepreg (min3)",perSD_hi],
             q=m3[method=="sleepreg (min3)",MF_Q5vQ1], ql=m3[method=="sleepreg (min3)",Q5_lo], qh=m3[method=="sleepreg (min3)",Q5_hi])
), use.names=TRUE)

ord <- c("Original (hand, published)","Hand (reproduction)","GGIR","sleepreg (min 5 d)","sleepreg (min 3 d)")
long <- rbind(
  dd[,.(method, panel="per +1 SD",            MF=perSD, lo=pl, hi=ph)],
  dd[,.(method, panel="top vs bottom quintile",MF=q,    lo=ql, hi=qh)])
long[, method := factor(method, levels=rev(ord))]
long[, panel := factor(panel, levels=c("per +1 SD","top vs bottom quintile"))]

p1 <- ggplot(long, aes(MF, method)) +
  geom_vline(xintercept=1, linetype="dashed", color="grey50") +
  geom_errorbar(aes(xmin=lo, xmax=hi), orientation="y", width=.22) +
  geom_point(size=2.4) +
  facet_wrap(~panel, scales="free_x") +
  labs(x="Multiplication factor on BMI (exp(coef)), 95% confidence interval",
       y=NULL,
       title="SRI-BMI effect size by SRI calculator (identical input, fully adjusted)") +
  theme_minimal(base_size=11) +
  theme(plot.title=element_text(size=11), panel.grid.minor=element_blank())
ggsave(file.path(figdir,"forest_MF.pdf"), p1, width=8.4, height=3.0)

# ---- agreement scatter vs hand ----
all <- as.data.table(readRDS(file.path(inter,"sri_all_methods.rds")))[!is.na(SRI_orig)]
sc <- rbind(
  all[is.finite(SRI_hand)&is.finite(SRI_GGIR),     .(hand=SRI_hand, pkg=SRI_GGIR,     which="GGIR")],
  all[is.finite(SRI_hand)&is.finite(SRI_sleepreg), .(hand=SRI_hand, pkg=SRI_sleepreg, which="sleepreg")])
p2 <- ggplot(sc, aes(hand, pkg)) +
  geom_point(alpha=.12, size=.5) +
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed") +
  geom_smooth(method="lm", se=FALSE, color="blue", linewidth=.7, formula=y~x) +
  facet_wrap(~which) + coord_equal() +
  labs(x="Hand SRI", y="Package SRI",
       title="Package SRI vs hand SRI on identical sleep/wake input (red = identity)") +
  theme_minimal(base_size=11) + theme(plot.title=element_text(size=11))
ggsave(file.path(figdir,"agreement_scatter.pdf"), p2, width=8.4, height=4.2)
cat("figures written to", figdir, "\n")
