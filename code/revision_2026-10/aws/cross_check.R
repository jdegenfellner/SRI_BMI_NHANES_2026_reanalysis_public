# cross_check.R <mac_base> <cloud_base> <liste.tsv> : vergleicht die auf dem Mac und die in der Cloud (Linux) berechneten Reihen derselben Teilnehmer
# GGIR: Epochenserie (alle Spalten), SWaN: Minuten-/30-s-Vorhersagen. Ergebnis: Zahl identischer Teilnehmer, Anteil abweichender Epochen, maximale Abweichung der Zahlenspalten.
suppressMessages(library(data.table))
a <- commandArgs(TRUE); mac <- a[1]; cl <- a[2]; L <- fread(a[3])
cmp <- function(f1, f2) { if (!file.exists(f1) || !file.exists(f2)) return(list(status = "Datei fehlt", n = NA, ndiff = NA, maxabs = NA))
  x <- fread(f1); y <- fread(f2); if (!identical(dim(x), dim(y)) || !identical(names(x), names(y))) return(list(status = "Dimension/Spalten verschieden", n = nrow(x), ndiff = NA, maxabs = NA))
  bad <- rep(FALSE, nrow(x)); mx <- 0
  for (j in names(x)) { u <- x[[j]]; v <- y[[j]]
    if (is.numeric(u) && is.numeric(v)) { d <- abs(u - v); d[is.na(u) & is.na(v)] <- 0; d[is.na(d)] <- Inf; mx <- max(mx, d[is.finite(d)], 0); bad <- bad | d > 1e-6 } else { e <- (is.na(u) & is.na(v)) | (!is.na(u) & !is.na(v) & as.character(u) == as.character(v)); bad <- bad | !e } }
  list(status = if (!any(bad)) "identisch" else "abweichend", n = nrow(x), ndiff = sum(bad), maxabs = mx) }
res <- rbindlist(lapply(seq_len(nrow(L)), function(i) { s <- L$seqn[i]; c <- L$cycle[i]
  g <- cmp(file.path(mac, "_ggir", c, s, paste0(s, "_epochs.csv.gz")), file.path(cl, "_ggir", c, s, paste0(s, "_epochs.csv.gz")))
  w <- cmp(file.path(mac, "_swan", c, paste0(s, "_final.csv.gz")), file.path(cl, "_swan", c, paste0(s, "_final.csv.gz")))
  data.table(seqn = s, cycle = c, ggir = g$status, ggir_n = g$n, ggir_ndiff = g$ndiff, ggir_maxabs = g$maxabs, swan = w$status, swan_n = w$n, swan_ndiff = w$ndiff) }))
print(res); fwrite(res, "cross_check_result.csv")
cat(sprintf("\nGGIR identisch: %d/%d, SWaN identisch: %d/%d; abweichende GGIR-Epochen gesamt %s von %s, SWaN %s von %s\n",
  sum(res$ggir == "identisch"), nrow(res), sum(res$swan == "identisch"), nrow(res), sum(res$ggir_ndiff, na.rm = TRUE), sum(res$ggir_n, na.rm = TRUE), sum(res$swan_ndiff, na.rm = TRUE), sum(res$swan_n, na.rm = TRUE)))
