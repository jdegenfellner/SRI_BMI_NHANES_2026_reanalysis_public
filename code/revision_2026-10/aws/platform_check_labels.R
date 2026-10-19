# Label-level comparison of the derived series of participants processed both on the reference laptop and on another machine.
# A participant whose laptop copy was itself transferred from that machine (cloud_ingested.log, no later laptop run) is a copy, not a
# duplicate, and is skipped. For every genuine duplicate: differing GGIR epochs in the label columns (class_id, SleepPeriodTime,
# invalidepoch), the largest differences in ACC (mg) and angle (degrees), differing SWaN predictions and the largest difference in the
# smoothed class probabilities. Usage: Rscript platform_check_labels.R <mac_base> <cloud_root> <ingest_log> <out_dir>
suppressMessages(library(data.table))
a <- commandArgs(TRUE); mac <- a[1]; root <- a[2]; ing <- fread(a[3], header = FALSE, col.names = c("time", "seqn", "cycle", "source")); outdir <- a[4]
ing[, t := as.POSIXct(time, format = "%Y-%m-%dT%H:%M:%S", tz = "Europe/Zurich")]
res <- list()
for (m in list.dirs(root, recursive = FALSE)) for (sub in c("out_check", "out_cc2", "prod", "validation", "out")) { d <- file.path(m, sub); if (!dir.exists(file.path(d, "_ggir"))) next
  for (cyc in c("pax_g", "pax_h")) { ids <- list.dirs(file.path(d, "_ggir", cyc), recursive = FALSE, full.names = FALSE); ids <- ids[grepl("^[0-9]+$", ids)]
    for (s in ids) { fm <- file.path(mac, "_ggir", cyc, s, paste0(s, "_epochs.csv.gz")); done <- file.path(mac, "_ggir", cyc, "done", paste0(s, ".txt")); if (!file.exists(fm) || !file.exists(done)) next
      src <- paste(basename(m), sub, sep = "/"); e <- ing[seqn == as.integer(s) & source == src]
      copied <- nrow(e) > 0 && as.numeric(difftime(file.mtime(done), max(e$t), units = "secs")) < 120   # laptop copy is the transferred cloud output
      if (copied) next
      x <- fread(fm); y <- fread(file.path(d, "_ggir", cyc, s, paste0(s, "_epochs.csv.gz")))
      if (!identical(dim(x), dim(y))) { res[[length(res) + 1]] <- data.table(machine = basename(m), folder = sub, seqn = s, cycle = cyc, note = "GGIR dimensions differ"); next }
      lab <- sum(x$class_id != y$class_id | x$SleepPeriodTime != y$SleepPeriodTime | x$invalidepoch != y$invalidepoch)
      acc <- max(abs(x$ACC - y$ACC), na.rm = TRUE); ang <- max(abs(x$angle - y$angle), na.rm = TRUE)
      w1 <- file.path(mac, "_swan", cyc, paste0(s, "_final.csv.gz")); w2 <- file.path(d, "_swan", cyc, paste0(s, "_final.csv.gz"))
      if (file.exists(w1) && file.exists(w2)) { u <- fread(w1); v <- fread(w2)
        if (identical(dim(u), dim(v))) { sw <- sum(u$PREDICTION != v$PREDICTION); pr <- max(abs(as.matrix(u[, .(PROB_WEAR_SMOOTH, PROB_SLEEP_SMOOTH, PROB_NWEAR_SMOOTH)]) - as.matrix(v[, .(PROB_WEAR_SMOOTH, PROB_SLEEP_SMOOTH, PROB_NWEAR_SMOOTH)])), na.rm = TRUE) } else { sw <- NA; pr <- NA } } else { sw <- NA; pr <- NA }
      res[[length(res) + 1]] <- data.table(machine = basename(m), folder = sub, seqn = s, cycle = cyc, note = "", ggir_epochs = nrow(x), ggir_label_diff = lab, acc_max_mg = acc, angle_max_deg = ang, swan_windows = if (is.na(sw)) NA_integer_ else nrow(u), swan_class_diff = sw, swan_prob_max = pr)
      cat(basename(m), sub, s, "labels", lab, "acc", acc, "angle", signif(ang, 3), "swan", sw, signif(pr, 3), "\n") } } }
out <- rbindlist(res, fill = TRUE); fwrite(out, file.path(outdir, "platform_check_labels.csv"))
cat("\n== genuine duplicates by machine and folder ==\n")
print(out[note == "", .(duplicates = .N, ggir_labels_identical = sum(ggir_label_diff == 0), max_label_diff = max(ggir_label_diff), acc_max_mg = max(acc_max_mg), angle_max_deg = signif(max(angle_max_deg), 3), swan_classes_identical = sum(swan_class_diff == 0, na.rm = TRUE), swan_prob_max = signif(max(swan_prob_max, na.rm = TRUE), 3)), by = .(machine, folder)][order(machine, folder)])
