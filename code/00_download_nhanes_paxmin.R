# =====================================================================
# 00_download_nhanes_paxmin.R
# Obtain the raw NHANES minute files PAXMIN_G (2011-12) and PAXMIN_H
# (2013-14) and convert them to the compact .RDS the pipeline (10_*.R)
# consumes. These files are ~1 GB each and cannot be shipped on GitHub;
# they are public-domain NHANES public-use data (see data/SOURCES.md).
#
# The small header files PAXHD_G/H are SHIPPED in the repo (~0.4 MB each),
# so only the two PAXMIN files need fetching.
#
# DOWNLOAD ROUTE: the nhanesA package works and is the recommended route -- this
# script tries it first (verified: nhanesA reaches wwwn.cdc.gov and returns the
# data). A plain curl / download.file is BLOCKED by CDC bot protection (it returns
# a ~20 KB HTML page instead of the file), so that fallback usually fails. If
# nhanesA is unavailable, download PAXMIN_G.XPT / PAXMIN_H.XPT manually from the
# documentation pages below into <root>/01_recompute_SRI/data_local/ and re-run;
# the script then converts each .XPT to the recoded .RDS.
#   install.packages("nhanesA")   # if needed
# =====================================================================
suppressMessages({library(haven)})
options(timeout = 60 * 60)

base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())
ddir <- file.path(base, "01_recompute_SRI", "data_local")   # matches 10_*.R
dir.create(ddir, recursive = TRUE, showWarnings = FALSE)
UA <- paste("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7)",
            "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120 Safari/537.36")

files <- list(
  list(seg="2011-2012", name="PAXMIN_G", rds="PAXMIN_2011_12.RDS",
       doc="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/PAXMIN_G.htm"),
  list(seg="2013-2014", name="PAXMIN_H", rds="PAXMIN_2013_14.RDS",
       doc="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/PAXMIN_H.htm")
)

# PAXPREDM recoding (integers 1-4 -> labels), matching the NCHS codebook:
# 1 = Wake wear, 2 = Sleep wear, 3 = Non wear, 4 = Unknown. See SOURCES.md.
recode_paxpredm <- function(x)
  factor(as.character(x), levels = c("1","2","3","4"),
         labels = c("wake","sleep","non-wear","unknown"))

get_raw <- function(f, xpt_path) {
  # (a) local .XPT already downloaded (e.g. manually) -> use it
  if (file.exists(xpt_path) && file.size(xpt_path) > 1e7)
    return(as.data.frame(haven::read_xpt(xpt_path)))
  # (b) try nhanesA (may work in some environments)
  if (requireNamespace("nhanesA", quietly = TRUE)) {
    message("Trying nhanesA for ", f$name, " ...")
    d <- try(nhanesA::nhanes(f$name, translated = FALSE), silent = TRUE)
    if (!inherits(d, "try-error") && is.data.frame(d) && nrow(d) > 1e5) return(d)
  }
  # (c) try a direct download (usually blocked by CDC bot protection)
  url <- sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s.XPT", f$seg, f$name)
  message("Trying direct download ", url, " ...")
  ok <- tryCatch(utils::download.file(url, xpt_path, mode="wb", quiet=FALSE,
                   headers=c("User-Agent"=UA)), error=function(e) 1L)
  if (ok == 0 && file.exists(xpt_path) && file.size(xpt_path) > 1e7)
    return(as.data.frame(haven::read_xpt(xpt_path)))
  NULL
}

for (f in files) {
  out <- file.path(ddir, f$rds)
  if (file.exists(out)) { message(f$rds, " already present, skipping."); next }
  dat <- get_raw(f, file.path(ddir, paste0(f$name, ".XPT")))
  if (is.null(dat)) {
    stop("Could not obtain ", f$name, " automatically (CDC bot protection).\n",
         "  Download it manually from:\n    ", f$doc, "\n",
         "  save ", f$name, ".XPT to:  ", ddir, "\n  then re-run this script.",
         call. = FALSE)
  }
  if ("PAXPREDM" %in% names(dat)) dat$PAXPREDM <- recode_paxpredm(dat$PAXPREDM)
  saveRDS(dat, out); message("  wrote ", f$rds, "  (", nrow(dat), " rows)")
  xpt <- file.path(ddir, paste0(f$name, ".XPT"))
  if (file.exists(xpt)) file.remove(xpt)   # drop the ~1 GB .XPT
}
message("\nDone. Next: run 10_build_common_series.R")
