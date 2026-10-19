suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())   # working tree root (contains 00_original_paper_materials/ and 01_recompute_SRI/)
an <- as.data.table(readRDS(file.path(base,"00_original_paper_materials/data/Combined_data_sets_2_waves_analytic_data_set.RDS")))[, .(SEQN, SRI_orig = SRI_no_imputation)]
s <- as.data.table(readRDS(file.path(base,"01_recompute_SRI/intermediate/series_1314.rds")))[SEQN %in% an$SEQN]
setorder(s, SEQN, UNIX_TIMESTAMP); s[, dayi := as.integer(PAXDAYM)]
# pick 3 participants with 7 consecutive full days at the 25th, 50th, 75th percentile of published SRI
nd <- s[, .(n = uniqueN(dayi), consec = all(diff(sort(unique(dayi))) == 1)), by = SEQN][n == 7 & consec]
cand <- merge(nd, an, by = "SEQN")[order(SRI_orig)]
pick <- cand[round(c(.25,.5,.75) * nrow(cand)), SEQN]
x <- s[SEQN %in% pick]
x[, day_start := UNIX_TIMESTAMP[1], by = .(SEQN, dayi)]; x[, tod_h := (UNIX_TIMESTAMP - day_start)/3600]
this <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE); here <- if (is.na(this)) getwd() else dirname(normalizePath(this))
src <- readLines(file.path(here, "variants4.R")); s0 <- grep("pass-2 operators", src); s1 <- grep("hand SRI \\(replicates", src); eval(parse(text = src[s0:(s1-1)]))
x[, bridged := v_bridgeNA(binary_sleep, 30L, 30L), by = SEQN]
x[, maj := v_thresh(binary_sleep, 60L, 0.5), by = SEQN]
x[, dayrel := dayi - min(dayi) + 1, by = SEQN]
lab <- merge(data.table(SEQN = pick), an, by = "SEQN")[match(pick, SEQN)]
x[, panel := factor(SEQN, levels = pick, labels = sprintf("Participant %s  (published SRI %.0f)", LETTERS[1:3], lab$SRI_orig))]
long <- melt(x[, .(panel, dayrel, tod_h, raw = binary_sleep, bridged, maj)], id.vars = c("panel","dayrel","tod_h"), variable.name = "series", value.name = "state")
long[, series := factor(series, levels = c("raw","bridged","maj"), labels = c("PAXPREDM as released\n(published analysis)", "Wake bouts <30 min bridged,\nisolated sleep bouts <30 min removed", "1-h rolling majority\n(\u226550% sleep \u2192 sleep)"))]
long[, st := factor(fifelse(is.na(state), "non-wear / unknown", fifelse(state == 1, "sleep", "wake")), levels = c("sleep","wake","non-wear / unknown"))]
p <- ggplot(long, aes(x = tod_h, y = dayrel, fill = st)) + geom_tile(height = 0.85, width = 1/60) +
  facet_grid(panel ~ series) + scale_y_reverse(breaks = 1:7, name = "Recording day") +
  scale_x_continuous(breaks = seq(0, 24, 6), labels = sprintf("%02d:00", seq(0,24,6)), name = "Clock time", expand = c(0,0)) +
  scale_fill_manual(values = c(sleep = "#1b3a6b", wake = "#f2e9d8", `non-wear / unknown` = "#c0c0c0"), name = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom", panel.grid = element_blank(), strip.text = element_text(size = 9))
ggsave("fig_examples_raw_vs_bridged.pdf", p, width = 12, height = 7.5); ggsave("fig_examples_raw_vs_bridged.png", p, width = 12, height = 7.5, dpi = 120)
# per-example stats
st <- x[, .(trans_raw = sum(diff(na.omit(binary_sleep)) != 0)/7, trans_br = sum(diff(na.omit(bridged)) != 0)/7, trans_maj = sum(diff(na.omit(maj)) != 0)/7,
            sleep_raw = sum(binary_sleep == 1, na.rm=TRUE)/7, sleep_br = sum(bridged == 1, na.rm=TRUE)/7, sleep_maj = sum(maj == 1, na.rm=TRUE)/7), by = panel]
print(st); cat("picked SEQN:", pick, "\nDONE\n")
