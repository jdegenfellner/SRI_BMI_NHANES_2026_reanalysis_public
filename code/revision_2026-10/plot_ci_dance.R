# Plot of the simulated confidence intervals (sim_ci_dance.R): one row per simulated sample, 95% (thick) and 99% (thin) bootstrap
# percentile intervals of the multiplication factor, point estimate as a dot, true factor as a solid line, 1 as a dashed grey line.
# Env: SRI_PROJECT_ROOT; reads sim_M2.csv and sim_M3w.csv from the sim_ci_dance results folder and writes ci_dance.pdf/.png there.
suppressMessages(library(data.table))
base <- Sys.getenv("SRI_PROJECT_ROOT"); dir <- file.path(base, "01_recompute_SRI/intermediate/revision_2026-10_results/sim_ci_dance")
panels <- list(M2 = "Model 2, fully adjusted, all participants", M3w = "Model 3, women")
draw <- function() {
  layout(matrix(c(1, 2, 3, 3), 2, byrow = TRUE), heights = c(10, 1.6)); par(mar = c(4.5, 3, 4.2, 1))
  for (m in names(panels)) {
    x <- fread(file.path(dir, paste0("sim_", m, ".csv"))); truth <- x$truth[1]; ns <- nrow(x)
    cov95 <- x[, lo95 <= truth & truth <= hi95]; cov99 <- x[, lo99 <= truth & truth <= hi99]
    xl <- range(c(x$lo99, x$hi99, truth, 1)); xl <- xl + c(-1, 1) * 0.06 * diff(xl)
    plot(NA, xlim = xl, ylim = c(0.5, ns + 0.5), yaxt = "n", ylab = "", xlab = "Multiplication factor on BMI, top versus bottom SRI quintile", main = "", bty = "l")
    axis(2, at = seq_len(ns), labels = seq_len(ns), las = 1, cex.axis = 0.75); mtext("simulated sample", side = 2, line = 2, cex = 0.85)
    abline(v = 1, lty = 2, col = "grey45"); abline(v = truth, lty = 1, col = "black")
    for (i in seq_len(ns)) {
      segments(x$lo99[i], i, x$hi99[i], i, lwd = 1, col = "grey55")
      segments(x$lo95[i], i, x$hi95[i], i, lwd = 3.2, col = if (cov95[i]) "#1f4e9c" else "#c0392b")
      points(x$mf[i], i, pch = 16, cex = 0.9, col = "black") }
    title(main = panels[[m]], cex.main = 1.05, line = 2.4)
    mtext(sprintf("true factor %.3f (solid line); point estimates %.3f to %.3f", truth, min(x$mf), max(x$mf)), side = 3, line = 1.3, cex = 0.75)
    mtext(sprintf("95%% intervals cover the true factor in %d of %d samples, 99%% intervals in %d of %d", sum(cov95), ns, sum(cov99), ns), side = 3, line = 0.4, cex = 0.75)
    cat(sprintf("%s: truth %.4f | point estimates %.3f-%.3f (SD %.4f) | mean width 95%% %.3f, 99%% %.3f | coverage 95%%: %d/%d, 99%%: %d/%d\n",
                m, truth, min(x$mf), max(x$mf), sd(x$mf), mean(x$hi95 - x$lo95), mean(x$hi99 - x$lo99), sum(cov95), ns, sum(cov99), ns))
  }
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", ncol = 2, legend = c("95% bootstrap interval (blue: contains the true factor; red: does not)", "point estimate", "99% bootstrap interval", ""),
         lwd = c(3.2, NA, 1, NA), pch = c(NA, 16, NA, NA), col = c("#1f4e9c", "black", "grey55", NA), bty = "n", cex = 0.85)
}
pdf(file.path(dir, "ci_dance.pdf"), width = 9.5, height = 6.2, pointsize = 12); draw(); invisible(dev.off())
png(file.path(dir, "ci_dance.png"), width = 9.5, height = 6.2, units = "in", res = 150, pointsize = 12); draw(); invisible(dev.off())
fig <- Sys.getenv("FIG_DIR", ""); if (nzchar(fig)) { pdf(file.path(fig, "fig_ci_dance.pdf"), width = 9.5, height = 6.2, pointsize = 12); draw(); invisible(dev.off()) }
cat("written:", file.path(dir, "ci_dance.pdf"), "\n")
