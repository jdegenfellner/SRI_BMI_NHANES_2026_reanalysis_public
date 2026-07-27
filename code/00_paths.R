# =====================================================================
# 00_paths.R -- resolve data & figure locations across the two layouts:
#   * public repo :  data/ , manuscript/figures/ , code/
#   * working tree:  01_recompute_SRI/... , 00_original_paper_materials/data/
# `base` must already be defined (Sys.getenv("SRI_PROJECT_ROOT", getwd())).
# Source this AFTER a script's own path definitions so it overrides them.
# =====================================================================
if (!exists("base")) base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())

.public <- dir.exists(file.path(base, "code")) &&
           dir.exists(file.path(base, "manuscript"))

if (.public) {
  rdir     <- base
  inter    <- file.path(base, "data")                # shipped derived data + outputs
  origdata <- file.path(base, "data")                # analytic modeling frames
  ddir     <- file.path(base, "data", "data_local")  # raw XPT/RDS (download target)
  figdir   <- file.path(base, "manuscript", "figures")
  sdir     <- file.path(base, "manuscript", "figures_suppl")
  src      <- file.path(base, "code", "sleepreg_src")
} else {
  rdir     <- file.path(base, "01_recompute_SRI")
  inter    <- file.path(rdir, "intermediate")
  origdata <- file.path(base, "00_original_paper_materials", "data")
  ddir     <- file.path(rdir, "data_local")
  figdir   <- file.path(base, "02_commentary_draft", "latex", "figures")
  sdir     <- file.path(base, "02_commentary_draft", "latex", "figures_suppl")
  src      <- file.path(rdir, "sleepreg_src")
}
for (.d in c(inter, figdir, sdir)) dir.create(.d, recursive = TRUE, showWarnings = FALSE)

# safety-net resolver: find an input file by basename across known data dirs
insrc <- function(name) {
  for (.d in unique(c(inter, origdata, ddir, base))) {
    f <- file.path(.d, name); if (file.exists(f)) return(f)
  }
  stop("Data file not found in any known location: ", name, call. = FALSE)
}
