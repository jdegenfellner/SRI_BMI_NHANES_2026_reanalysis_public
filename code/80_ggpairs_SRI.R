# =====================================================================
# Three-way (and four-way) SRI agreement: GGally::ggpairs with all-vs-all
# point clouds. Upper panels show PEARSON r with 95% CI (cor.test) --
# NO significance stars.
# =====================================================================
suppressMessages({library(data.table); library(ggplot2); library(GGally)})
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())
source(file.path(base, if (dir.exists(file.path(base, "code"))) "code" else "01_recompute_SRI", "00_paths.R"))
# figdir provided by 00_paths.R

all <- as.data.table(readRDS(file.path(inter,"sri_all_methods.rds")))[!is.na(SRI_orig)]

# upper panel: Pearson r with 95% CI, no stars
upper_rci <- function(data, mapping, ...){
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  ok <- is.finite(x) & is.finite(y)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "pearson"))
  lab <- sprintf("r = %.2f\n[%.2f, %.2f]\nn = %d",
                 ct$estimate, ct$conf.int[1], ct$conf.int[2], sum(ok))
  ggplot(data.frame(x=1,y=1), aes(x,y)) +
    annotate("text", x=1, y=1, label=lab, size=4.3, lineheight=0.95) +
    theme_void()
}
lower_pts <- function(data, mapping, ...){
  ggplot(data, mapping) +
    geom_point(alpha=0.10, size=0.4) +
    geom_abline(slope=1, intercept=0, colour="red", linetype="dashed", linewidth=0.4) +
    theme_minimal(base_size=10)
}
diag_dens <- function(data, mapping, ...){
  ggplot(data, mapping) + geom_density(fill="grey80", colour="grey30", linewidth=0.4) +
    theme_minimal(base_size=10)
}

make <- function(cols, labs, file, w, h){
  df <- as.data.frame(all[, ..cols]); colnames(df) <- labs
  p <- ggpairs(df,
               lower = list(continuous = lower_pts),
               diag  = list(continuous = diag_dens),
               upper = list(continuous = upper_rci),
               title = "Sleep Regularity Index across calculators (identical input): point clouds, Pearson r with 95% CI") +
    theme(plot.title = element_text(size = 11))
  ggsave(file.path(figdir, file), p, width=w, height=h)
  cat("wrote", file, "\n")
}

# three calculators (Published, GGIR, sleepreg) -- used as Figure 1 in the main text
make(c("SRI_orig","SRI_GGIR","SRI_sleepreg"),
     c("Published","GGIR","sleepreg"),
     "ggpairs_SRI_3way.pdf", 8.4, 8.0)
cat("DONE\n")
