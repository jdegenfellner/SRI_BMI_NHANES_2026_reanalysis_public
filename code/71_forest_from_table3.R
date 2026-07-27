# Forest plot consistent with the Table-3-comparable results (Model 2 full adj.
# and Model 3 women), one row per SRI calculator, with bootstrap 95% CIs.
suppressMessages({library(data.table); library(ggplot2)})
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())
source(file.path(base, if (dir.exists(file.path(base, "code"))) "code" else "01_recompute_SRI", "00_paths.R"))
# figdir provided by 00_paths.R
t <- as.data.table(readRDS(file.path(inter,"table3_replication.rds")))
t <- t[calculator != "hand"]   # baseline = published SRI ("orig"), no hand re-implementation

# parse "0.950 (0.931; 0.969)"  and  "Female: 0.966 (0.937; 0.994)"
p3 <- function(s){ m <- regmatches(s, regexec("([0-9.]+) \\(([0-9.]+); ([0-9.]+)\\)", s))[[1]]
                   as.numeric(m[2:4]) }
lab <- c(orig="Published", GGIR="GGIR",
         sleepreg5="sleepreg (min 5 d)", sleepreg3="sleepreg (min 3 d)")
rows <- rbindlist(lapply(1:nrow(t), function(i){
  m2 <- p3(t$M2_fully_adj[i])
  fem <- p3(t$M3_sex_other[i])   # M3_sex_other = Female (ref level is Male)
  rbind(
    data.table(method=lab[t$calculator[i]], panel="Fully adjusted (Model 2)", MF=m2[1], lo=m2[2], hi=m2[3]),
    data.table(method=lab[t$calculator[i]], panel="Women, fully adjusted (Model 3)", MF=fem[1], lo=fem[2], hi=fem[3]))
}))
ord <- rev(c("Published","GGIR","sleepreg (min 5 d)","sleepreg (min 3 d)"))
rows[, method := factor(method, levels=ord)]
rows[, panel := factor(panel, levels=c("Fully adjusted (Model 2)","Women, fully adjusted (Model 3)"))]

p <- ggplot(rows, aes(MF, method)) +
  geom_vline(xintercept=1, linetype="dashed", colour="grey50") +
  geom_errorbar(aes(xmin=lo, xmax=hi), orientation="y", width=.22) +
  geom_point(size=2.4) +
  facet_wrap(~panel) +
  labs(x="Multiplication factor on BMI, Q1->Q5 (95% confidence interval)", y=NULL,
       title="SRI-BMI effect size by calculator -- published Table-3 method (identical models)") +
  theme_minimal(base_size=11) + theme(plot.title=element_text(size=10.5), panel.grid.minor=element_blank())
ggsave(file.path(figdir,"forest_MF.pdf"), p, width=8.6, height=2.9)
cat("forest_MF.pdf rebuilt from Table-3 results\n")
