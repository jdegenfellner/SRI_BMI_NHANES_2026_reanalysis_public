# Compact backup of the derived minute-level sleep/wake series (GGIR SIB, GGIR SPT, SWaN) of every finished participant.
# One row per participant and platform; each series is a string with one character per minute (1 sleep, 0 wake, . missing),
# minute 0 = start of the raw recording. The raw GGIR/SWaN output stays on the external disk; this file is the small, versionable result.
#
# Environment variables (all optional):
#   RAW_ROOT   local run (platform "mac"), default /Volumes/LaCie/NHANES_1112_1314
#   CLOUD_DIR  staging folder of the cloud instances (<instance>/<run>/_ggir, _swan), default <RAW_ROOT>/_cloud
#   OUT_FILE   default data/raw_derivation/minute_series.tsv.gz (relative to the working directory = repository root)
# Rows already present in OUT_FILE are kept; only new (SEQN, platform) pairs are computed.
suppressMessages(library(data.table))
raw0   <- Sys.getenv("RAW_ROOT", "/Volumes/LaCie/NHANES_1112_1314")
cloud  <- Sys.getenv("CLOUD_DIR", file.path(raw0, "_cloud"))
outf   <- Sys.getenv("OUT_FILE", "data/raw_derivation/minute_series.tsv.gz")
dir.create(dirname(outf), showWarnings = FALSE, recursive = TRUE)
an <- unique(as.data.table(readRDS("data/df_an1_imputed_MAIN.RDS"))$SEQN)
this <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE)   # Rscript encodes spaces in --file as ~+~
here <- if (is.na(this)) getwd() else dirname(normalizePath(this))
src <- readLines(file.path(here, "raw_derivation_participants.R"))
a <- grep("^read_meta <- ", src); b <- grep("^# fragmentation on the retained days", src) - 1
raw <- raw0; eval(parse(text = src[a:b]))                       # read_meta, lab2bin, ggir_minutes, swan_minutes (identical code to the analysis)

cl <- unlist(lapply(list.dirs(cloud, recursive = FALSE), function(d) { s <- list.dirs(d, recursive = FALSE); s <- s[dir.exists(file.path(s, "_ggir"))]; setNames(s, paste(basename(d), basename(s), sep = "/")) }))
roots <- c(list(mac = raw0), as.list(cl))
old <- if (file.exists(outf)) fread(outf, colClasses = list(character = c("sib", "spt", "swan"))) else data.table()
have <- if (nrow(old)) paste(old$SEQN, old$platform) else character()
ids <- function(r, tool, cyc) as.integer(sub("\\.txt$", "", list.files(file.path(r, paste0("_", tool), cyc, "done"), pattern = "\\.txt$")))
enc <- function(dt, col, L) { v <- rep(".", L); j <- dt$minute + 1L; x <- dt[[col]]; ok <- j >= 1L & j <= L & !is.na(x); v[j[ok]] <- as.character(x[ok]); paste(v, collapse = "") }
ing <- file.path(raw0, "_work", "cloud_ingested.log")
ingested <- if (file.exists(ing)) as.integer(sapply(strsplit(readLines(ing), " "), `[`, 2)) else integer()   # already exported from the cloud staging folder
new <- list()
for (pf in names(roots)) {
  r <- roots[[pf]]; raw <- r
  for (cyc in c("pax_g", "pax_h")) {
    both <- intersect(ids(r, "ggir", cyc), ids(r, "swan", cyc)); both <- both[both %in% an]
    for (sq in both) {
      if (paste(sq, pf) %in% have || (pf == "mac" && sq %in% ingested)) next
      row <- tryCatch({
        gm <- ggir_minutes(sq, cyc); sm <- swan_minutes(sq, cyc); if (is.null(gm) || is.null(sm)) NULL else {
          L <- max(gm$minute, sm$minute) + 1L
          data.table(SEQN = sq, cycle = cyc, platform = pf, n_min = L, sib = enc(gm, "GGIR_SIB", L), spt = enc(gm, "GGIR_SPT", L), swan = enc(sm, "SWaN", L)) } },
        error = function(e) { message("skip ", sq, " ", pf, ": ", conditionMessage(e)); NULL })
      if (!is.null(row)) new[[length(new) + 1]] <- row
    }
  }
}
res <- rbindlist(c(list(old), new), fill = TRUE); setorder(res, SEQN, platform)
tmp <- paste0(outf, ".tmp"); fwrite(res, tmp, sep = "\t", compress = "gzip"); file.rename(tmp, outf)
cat(sprintf("%s: %d new rows, %d rows in total (%d participants), file %.1f MB\n", format(Sys.time(), "%FT%T"), length(new), nrow(res), uniqueN(res$SEQN), file.size(outf) / 1e6))
print(res[, .N, by = platform])
