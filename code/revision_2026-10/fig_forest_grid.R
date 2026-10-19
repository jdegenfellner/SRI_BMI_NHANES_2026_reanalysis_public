suppressMessages({library(data.table); library(ggplot2)})
t <- readRDS("na_models_all.rds")[calculator %in% c("SRI_hand","SRI_GGIR","SRI_sleepregE5","SRI_sleepregE3")]
pr <- function(x) { m <- regmatches(x, regexec("([0-9.]+) \\(([0-9.]+), ([0-9.]+)\\)", x)); t(sapply(m, function(z) as.numeric(z[2:4]))) }
a <- pr(t$M2); w <- pr(sub("^Female ", "", t$M3_other))
d <- rbind(data.table(t[, .(variant, calculator)], model = "Model 2, fully adjusted", est = a[,1], lo = a[,2], hi = a[,3]),
           data.table(t[, .(variant, calculator)], model = "Model 3, women", est = w[,1], lo = w[,2], hi = w[,3]))
d[, variant := factor(variant, levels = rev(c("raw","maj1h50","ggir1h80","bridge2","spt2")), labels = rev(c("Released classification, unsmoothed","Rolling majority, 1 h","GGIR SRI1 smoothing, 1 h, 80%","Awakenings <30 min bridged","Main sleep period only")))]
d[, calculator := factor(calculator, levels = c("SRI_hand","SRI_GGIR","SRI_sleepregE5","SRI_sleepregE3"), labels = c("Phillips formula","GGIR","sleepreg, 5 days","sleepreg, 3 days"))]
pub <- data.table(model = c("Model 2, fully adjusted","Model 3, women"), est = c(0.950, 0.918))
p <- ggplot(d, aes(x = est, y = variant, colour = calculator, shape = calculator)) + geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_vline(data = pub, aes(xintercept = est), linetype = 3, colour = "#c0392b") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.7), size = 0.35) + facet_wrap(~model) +
  scale_colour_manual(values = c("#1b3a6b","#2e86ab","#e07b39","#b5651d"), name = NULL) + scale_shape_manual(values = c(16,17,15,18), name = NULL) +
  labs(x = "Multiplication factor on BMI, top vs bottom SRI quintile (95% bootstrap CI)", y = NULL) + theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave("fig_forest_grid.pdf", p, width = 10, height = 5.2); cat("forest figure written\n")
