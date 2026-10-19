suppressMessages({library(data.table); library(ggplot2)})
byh <- readRDS("byh_true.rds"); fs <- readRDS("frag_summary.rds")
p1 <- ggplot(byh, aes(hour + 0.5, p_sleep)) + geom_col(fill = "#1b3a6b", width = 0.9) +
  scale_x_continuous(breaks = seq(0, 24, 3), labels = sprintf("%02d:00", seq(0, 24, 3)), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "Clock time", y = "Share of minutes classified as sleep", title = "A. PAXPREDM sleep share by clock hour (n = 7085, all retained days)") + theme_minimal(base_size = 11)
sb <- fs$sb; sb[, cl := cut(len, c(0, 5, 15, 30, 60, 120, 1e4), labels = c("1-5", "6-15", "16-30", "31-60", "61-120", ">120"))]
d2 <- sb[, .(bouts = .N, minutes = sum(len)), by = cl][, `:=`(share_bouts = bouts/sum(bouts), share_min = minutes/sum(minutes))]
d2l <- melt(d2[, .(cl, `share of sleep bouts` = share_bouts, `share of sleep minutes` = share_min)], id.vars = "cl")
p2 <- ggplot(d2l, aes(cl, value, fill = variable)) + geom_col(position = "dodge") + scale_fill_manual(values = c("#1b3a6b", "#9bb0d1"), name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) + labs(x = "Sleep bout length (min)", y = NULL, title = "B. Length of uninterrupted sleep bouts in PAXPREDM") + theme_minimal(base_size = 11) + theme(legend.position = "bottom")
library(patchwork); p <- p1 / p2
ggsave("fig_input_diagnostics.pdf", p, width = 8, height = 8); ggsave("fig_input_diagnostics.png", p, width = 8, height = 8, dpi = 110)
print(d2[, .(cl, bouts, share_bouts = round(share_bouts, 3), share_min = round(share_min, 3))]); cat("DONE\n")
