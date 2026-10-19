# Compares, for every cloud or laptop folder, the GGIR epoch series and SWaN predictions of the participants that were also processed
# on the main laptop (the reference), with the comparison of cross_check.R. Usage: Rscript platform_check_all.R <mac_base> <cloud_root> <out_dir>
#   mac_base   folder with _ggir/ and _swan/ of the main laptop;  cloud_root  folder with one subfolder per machine (out_check, out_cc2, prod, validation inside)
suppressMessages(library(data.table))
a <- commandArgs(TRUE); mac <- a[1]; root <- a[2]; outdir <- a[3]; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
cmp <- function(f1, f2) { if (!file.exists(f1) || !file.exists(f2)) return(list(status = "file missing", n = NA, ndiff = NA, maxabs = NA))
  x <- fread(f1); y <- fread(f2); if (!identical(dim(x), dim(y)) || !identical(names(x), names(y))) return(list(status = "dimensions differ", n = nrow(x), ndiff = NA, maxabs = NA))
  bad <- rep(FALSE, nrow(x)); mx <- 0
  for (j in names(x)) { u <- x[[j]]; v <- y[[j]]
    if (is.numeric(u) && is.numeric(v)) { d <- abs(u - v); d[is.na(u) & is.na(v)] <- 0; d[is.na(d)] <- Inf; mx <- max(mx, d[is.finite(d)], 0); bad <- bad | d > 1e-6 } else { e <- (is.na(u) & is.na(v)) | (!is.na(u) & !is.na(v) & as.character(u) == as.character(v)); bad <- bad | !e } }
  list(status = if (!any(bad)) "identical" else "different", n = nrow(x), ndiff = sum(bad), maxabs = mx) }
res <- list()
for (m in list.dirs(root, recursive = FALSE)) for (sub in c("out_check", "out_cc2", "prod", "validation", "out")) { d <- file.path(m, sub); if (!dir.exists(file.path(d, "_ggir"))) next
  for (cyc in c("pax_g", "pax_h")) { ids <- list.dirs(file.path(d, "_ggir", cyc), recursive = FALSE, full.names = FALSE); ids <- ids[grepl("^[0-9]+$", ids)]
    for (s in ids) { fm <- file.path(mac, "_ggir", cyc, s, paste0(s, "_epochs.csv.gz")); if (!file.exists(fm)) next   # only participants also processed on the reference laptop
      g <- cmp(fm, file.path(d, "_ggir", cyc, s, paste0(s, "_epochs.csv.gz")))
      w <- cmp(file.path(mac, "_swan", cyc, paste0(s, "_final.csv.gz")), file.path(d, "_swan", cyc, paste0(s, "_final.csv.gz")))
      res[[length(res) + 1]] <- data.table(machine = basename(m), folder = sub, seqn = s, cycle = cyc, ggir = g$status, ggir_n = g$n, ggir_ndiff = g$ndiff, ggir_maxabs = g$maxabs, swan = w$status, swan_n = w$n, swan_ndiff = w$ndiff)
      cat(basename(m), sub, s, g$status, w$status, "\n") } } }
out <- rbindlist(res); fwrite(out, file.path(outdir, "platform_check_all.csv")); print(out)
cat("\nby machine:\n"); print(out[, .(participants = .N, ggir_identical = sum(ggir == "identical"), swan_identical = sum(swan == "identical"), ggir_epochs_differing = sum(ggir_ndiff, na.rm = TRUE), max_abs = max(ggir_maxabs, na.rm = TRUE)), by = machine])
