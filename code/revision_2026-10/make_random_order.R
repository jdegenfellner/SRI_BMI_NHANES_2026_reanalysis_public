# Random processing order of the published sample for the raw-signal derivation (analytic_random_order.tsv).
# The 7085 SEQN of the analytic sample are sorted and permuted once with R's default generator and the seed 20261016;
# the rank is the position in this permutation. Env: SRI_PROJECT_ROOT (working-tree root with 00_original_paper_materials/).
suppressMessages(library(data.table))
this <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE)
here <- if (is.na(this)) getwd() else dirname(normalizePath(this))
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())
d <- as.data.table(readRDS(file.path(base, "00_original_paper_materials", "data", "df_an1_imputed_MAIN.RDS")))[, .(seqn = SEQN)]
d[, cycle := fifelse(seqn < 73557, "pax_g", "pax_h")]   # NHANES 2011-2012 SEQN 62161-71916, 2013-2014 SEQN 73557-83731
setorder(d, seqn); set.seed(20261016); ord <- d[sample(.N)]; ord[, rank := .I]
out <- ord[, .(rank, seqn, cycle)]
ref <- file.path(here, "analytic_random_order.tsv")
if (file.exists(ref)) { old <- fread(ref); cat("identical to the shipped file:", all(old$rank == out$rank) && all(old$seqn == out$seqn) && all(old$cycle == out$cycle), "\n") }
fwrite(out, file.path(Sys.getenv("OUT_DIR", getwd()), "analytic_random_order_regenerated.tsv"), sep = "\t")
