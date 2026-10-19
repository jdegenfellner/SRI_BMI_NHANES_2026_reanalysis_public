# Figure: distribution of the SRI (Phillips formula, common mask) for the released classification and the three series derived from the raw signal.
# Usage: Rscript fig_sri_by_derivation.R <raw_interim_folder> <output_folder>
suppressMessages({library(data.table); library(ggplot2)})
a <- commandArgs(TRUE); indir <- a[1]; outdir <- a[2]
res <- readRDS(file.path(indir, "raw_derivation_participants.rds"))
ord <- fread(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))), "analytic_random_order.tsv"))
rmax <- as.numeric(Sys.getenv("R_MAX", "1000")); samp <- ord[rank <= rmax, seqn]
r <- res[SEQN %in% samp]; full <- r[, .N, by = SEQN][N == 4, SEQN]; r <- r[SEQN %in% full]
labs <- c(PAXPREDM = "PAXPREDM as released", GGIR_SIB = "GGIR SIB (inactivity bouts in the sleep period)", GGIR_SPT = "GGIR SPT (sleep period)", SWaN = "SWaN rerun on the raw signal")
r[, variant := factor(labs[deriv], levels = labs)]
med <- r[, .(m = median(sri_A, na.rm = TRUE)), by = variant]
p <- ggplot(r, aes(sri_A)) + geom_histogram(binwidth = 2, fill = "#1b3a6b", colour = NA) + geom_vline(data = med, aes(xintercept = m), colour = "#c0392b") +
  geom_text(data = med, aes(x = m, y = Inf, label = sprintf("median %.1f", m)), vjust = 1.5, hjust = -0.05, colour = "#c0392b", size = 3.2) +
  facet_wrap(~variant, ncol = 1, scales = "free_y") + xlim(-30, 100) + labs(x = sprintf("SRI (Phillips formula, common minutes; n = %d participants)", length(full)), y = "Participants") + theme_minimal(base_size = 11)
ggsave(file.path(outdir, "fig_sri_by_derivation.pdf"), p, width = 7, height = 7.5); ggsave(file.path(outdir, "fig_sri_by_derivation.png"), p, width = 7, height = 7.5, dpi = 110); cat("DONE", length(full), "\n")
