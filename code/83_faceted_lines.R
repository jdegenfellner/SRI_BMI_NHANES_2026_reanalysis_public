# =====================================================================
# Cleaner version of the "several slopes" figure: facet by valid-days so the
# clouds don't overlap. Each panel: points + identity (dashed) + OLS fit
# (solid) + LOESS (dotted) + annotated slope & r. The poor linear fit at
# 5-6 days is the point -- agreement BREAKS DOWN, not just shifts.
# =====================================================================
suppressMessages({library(data.table); library(ggplot2)})
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())
source(file.path(base, if (dir.exists(file.path(base, "code"))) "code" else "01_recompute_SRI", "00_paths.R"))

all  <- as.data.table(readRDS(insrc("sri_all_methods.rds")))[!is.na(SRI_orig)]
feat <- as.data.table(read.csv(insrc("participant_ndays.csv")))[, .(SEQN, n_days)]
d <- merge(all[,.(SEQN,SRI_orig,SRI_GGIR,SRI_sleepreg)], feat, by="SEQN")[n_days %in% 5:7]

# long: one row per (participant, package)
L <- rbind(
  d[is.finite(SRI_sleepreg), .(SEQN, hand=SRI_orig, pkg=SRI_sleepreg, calc="sleepreg", n_days)],
  d[is.finite(SRI_GGIR),     .(SEQN, hand=SRI_orig, pkg=SRI_GGIR,     calc="GGIR",     n_days)])
L[, days := factor(n_days, levels=7:5, labels=c("7 days (consecutive)","6 days","5 days"))]
L[, calc := factor(calc, levels=c("sleepreg","GGIR"))]

# per-panel slope & r for annotation
ann <- L[, .(slope=coef(lm(pkg~hand))[2], r=cor(hand,pkg), n=.N), by=.(calc,days)]
ann[, lab := sprintf("slope = %.2f\nr = %.2f\nn = %d", slope, r, n)]

p <- ggplot(L, aes(hand, pkg)) +
  geom_point(alpha=.14, size=.45, colour="grey25") +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="red") +
  geom_smooth(method="loess", se=FALSE, colour="orange", linewidth=.7, linetype="dotted", formula=y~x) +
  geom_smooth(method="lm",    se=TRUE,  colour="blue",   linewidth=.8, formula=y~x) +
  geom_text(data=ann, aes(x=8, y=92, label=lab), hjust=0, vjust=1, size=3.1, colour="black") +
  facet_grid(calc ~ days) +
  coord_cartesian(xlim=c(0,100), ylim=c(-20,100)) +
  labs(x="Published SRI", y="Package SRI",
       title="Agreement with the published SRI by number of valid days (red = identity, blue = OLS, orange = LOESS)") +
  theme_minimal(base_size=10) + theme(plot.title=element_text(size=10))
ggsave(file.path(figdir,"faceted_lines_by_days.pdf"), p, width=10, height=6.4)
cat("=== per-panel fits ===\n"); print(ann[order(calc,days), .(calc,days,slope=round(slope,3),r=round(r,3),n)])
cat("figure written: faceted_lines_by_days.pdf\nDONE\n")
