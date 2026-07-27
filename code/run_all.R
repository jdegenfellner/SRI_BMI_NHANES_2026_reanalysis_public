#!/usr/bin/env Rscript
# =====================================================================
# run_all.R  --  reproduce every table and figure in the manuscript.
#
#   # default: reproduce everything from the shipped data/ (no download)
#   Rscript code/run_all.R
#
#   # full pipeline from raw NHANES (needs the ~2 GB PAXMIN download and the
#   # sleepreg source; see README) -- rebuilds the SRIs before Stage 2:
#   FROM_RAW=TRUE Rscript code/run_all.R
#
# Figures land in manuscript/figures and manuscript/figures_suppl; result
# tables/CSVs land in data/. Each script is self-contained and resolves the
# repo root from SRI_PROJECT_ROOT (set below).
# =====================================================================
t0 <- proc.time()[["elapsed"]]

# --- resolve repo root (folder that contains code/ and manuscript/) ---
.args <- commandArgs(trailingOnly = FALSE)
.fa   <- sub("^--file=", "", .args[grep("^--file=", .args)])
root  <- if (length(.fa)) { normalizePath(file.path(dirname(.fa[1]), "..")) } else { Sys.getenv("SRI_PROJECT_ROOT", getwd()) }
Sys.setenv(SRI_PROJECT_ROOT = root)
code  <- file.path(root, "code")
from_raw <- toupper(Sys.getenv("FROM_RAW", "FALSE")) %in% c("TRUE", "1", "YES")

# Each script runs in its OWN Rscript process: a clean environment and graphics
# device per step, exactly as the scripts are designed to run standalone. The
# child inherits SRI_PROJECT_ROOT (set above), so every script resolves the repo.
run <- function(script, note = "") {
  cat(sprintf("\n=== %s %s ===\n", script, note))
  status <- system2("Rscript", shQuote(file.path(code, script)), stdout = "", stderr = "")
  ok <- identical(status, 0L)
  if (!ok) cat(sprintf("  !! FAILED (exit %s)\n", status))
  ok
}

cat("Repo root :", root, "\nFROM_RAW  :", from_raw, "\n")

# --- Stage 1 (optional): rebuild the three SRIs from the raw NHANES series ---
if (from_raw) {
  cat("\n##### STAGE 1: rebuild SRIs from raw NHANES (download + compute) #####\n")
  for (s in c("00_download_nhanes_paxmin.R", "10_build_common_series.R",
              "11_participant_features.R", "20_compute_three_SRI.R",
              "60_sleepreg_minSRIdays3.R")) run(s, "[raw]")
}

# --- Stage 2 (always): models, tables and figures from the derived data ---
cat("\n##### STAGE 2: models, tables and figures (from shipped data/) #####\n")
core <- c(
  "01_published_sri.R",          # check the published-SRI snapshot == the analysis baseline
  "30_run_models.R",             # survey-weighted BMI models per calculator -> model_comparison
  "31_table3_replication.R",     # Table 2 (MF per calculator, published Table-3 method)
  "70_figures_commentary.R",     # coefficient / MF overview figures
  "71_forest_from_table3.R",     # Figure 2A  (forest of MFs)
  "74_forest_diff.R",            # Figure 2B  (paired-bootstrap MF differences)
  "76_figure2_by_calculator.R",  # Figure 3   (original Fig 2 replicated per calculator)
  "83_faceted_lines.R",          # Figure 4   (agreement by valid-days)
  "81_cloud_structure.R",        # Figure 5   (cloud structure of the disagreement)
  "72_subgroup_forest.R",        # Figure 6A  (MF by recording completeness)
  "75_subgroup_diff.R",          # Figure 6B  (paired-bootstrap subgroup differences)
  "80_ggpairs_SRI.R",            # Figure 1   (pairwise agreement)
  "40_agreement_figures.R",      # agreement summaries (Table 1 numbers, Bland-Altman)
  "96_supplement_plots.R",       # supplement: distribution, Bland-Altman, predicted BMI
  "97_model_assumptions.R",      # supplement: quantitative assumption table
  "98_check_model.R"             # supplement: check_model() panels + car::qqPlot()
)
res <- vapply(core, run, logical(1))

# verification table (synthetic input; needs the sleepreg source, which is LGPL
# and not bundled -- see README on how to point sleepreg_src at it)
has_sleepreg <- any(dir.exists(file.path(root, c("01_recompute_SRI/sleepreg_src", "code/sleepreg_src"))))
if (has_sleepreg) {
  res <- c(res, `95_verify_packages.R` = run("95_verify_packages.R", "[synthetic]"))
} else {
  cat("\n=== 95_verify_packages.R SKIPPED (sleepreg source not found; see README) ===\n")
}

cat(sprintf("\n##### DONE: %d/%d scripts ok in %.0fs #####\n",
            sum(res), length(res), proc.time()[["elapsed"]] - t0))
if (!all(res)) cat("FAILED:", paste(names(res)[!res], collapse = ", "), "\n")
