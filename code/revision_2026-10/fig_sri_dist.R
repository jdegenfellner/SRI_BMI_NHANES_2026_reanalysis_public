suppressMessages({library(data.table); library(ggplot2)})
r <- as.data.table(readRDS("variants4_sri.rds"))
vars <- c("raw","maj1h50","ggir1h80","bridge2","spt2")
labs <- c(raw = "Released classification, unsmoothed", maj1h50 = "Rolling majority, 1 h", ggir1h80 = "GGIR SRI1 smoothing, 1 h, 80%", bridge2 = "Awakenings <30 min bridged", spt2 = "Main sleep period only")
l <- melt(r[, c("SEQN", paste0("SRI_", vars)), with = FALSE], id.vars = "SEQN")[, variant := factor(sub("SRI_", "", variable), levels = vars, labels = labs[vars])]
med <- l[, .(m = median(value, na.rm = TRUE)), by = variant]
p <- ggplot(l, aes(value)) + geom_histogram(binwidth = 2, fill = "#1b3a6b", colour = NA) + geom_vline(data = med, aes(xintercept = m), colour = "#c0392b") +
  geom_text(data = med, aes(x = m, y = Inf, label = sprintf("median %.1f", m)), vjust = 1.5, hjust = -0.05, colour = "#c0392b", size = 3.2) +
  facet_wrap(~variant, ncol = 1, scales = "free_y") + xlim(-40, 100) + labs(x = "SRI (Phillips formula on each input variant)", y = "Participants") + theme_minimal(base_size = 11)
ggsave("fig_sri_by_input.pdf", p, width = 7, height = 9); ggsave("fig_sri_by_input.png", p, width = 7, height = 9, dpi = 110); cat("DONE\n")
