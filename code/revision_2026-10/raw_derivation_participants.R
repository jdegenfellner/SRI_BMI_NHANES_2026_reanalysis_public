# One row per participant and series: SRI, missingness, fragmentation and agreement with the released classification (PAXPREDM)
# for the sleep/wake series derived from the 80 Hz raw signal (GGIR SIB, GGIR SPT, SWaN).
#
# Environment variables (all optional):
#   SRI_PROJECT_ROOT  working-tree root (contains 00_original_paper_materials/ and 01_recompute_SRI/), default: working directory
#   RAW_ROOT          folder with _ggir/ and _swan/ outputs (default /Volumes/LaCie/NHANES_1112_1314)
#   PAXMIN_1112, PAXMIN_1314   PAXMIN data per cycle (RDS with SEQN, PAXDAYM, PAXSSNMP, PAXTSM, PAXQFM, PAXAISMM, PAXPREDM);
#                     default: 01_recompute_SRI/data_local/PAXMIN_2011_12.RDS and PAXMIN_2013_14.RDS
#   OUT_DIR           output folder (default: working directory)
#   MAX_N             process at most this many participants (default all)
#
# Series definitions (Methods, "Derivation of the sleep/wake series from the raw signal"):
#   A (main analysis)  common mask: only minutes that are wake or sleep in the rebuilt series, and wake or sleep in the derived series
#   B (sensitivity)    own missing minutes on the retained days of the rebuilt series; participants with >= 30% missing are excluded
# The minute index of the released data is PAXSSNMP / 4800 (sample counter of the 80 Hz signal); the lag for the SRI is 1440 minutes
# on this grid, so the SRI equals the Phillips formula on the continuous timeline with missing minutes as NA.
suppressMessages(library(data.table))
setDTthreads(2)
base   <- Sys.getenv("SRI_PROJECT_ROOT", getwd())
raw    <- Sys.getenv("RAW_ROOT", "/Volumes/LaCie/NHANES_1112_1314")
outdir <- Sys.getenv("OUT_DIR", getwd())
max_n  <- as.integer(Sys.getenv("MAX_N", "0"))
paxrds <- c(pax_g = Sys.getenv("PAXMIN_1112", file.path(base, "01_recompute_SRI/data_local/PAXMIN_2011_12.RDS")),
            pax_h = Sys.getenv("PAXMIN_1314", file.path(base, "01_recompute_SRI/data_local/PAXMIN_2013_14.RDS")))
t_start <- Sys.time()

an <- unique(as.data.table(readRDS(file.path(base, "00_original_paper_materials/data/df_an1_imputed_MAIN.RDS")))$SEQN)

# ---- which participants have which derived series ---------------------------------------------------------------------------------
done_ids <- function(tool, cyc) as.integer(sub("\\.txt$", "", list.files(file.path(raw, paste0("_", tool), cyc, "done"), pattern = "\\.txt$")))
avail <- rbindlist(lapply(c("pax_g", "pax_h"), function(cyc) rbindlist(list(
  data.table(SEQN = done_ids("ggir", cyc), cycle = cyc, tool = "ggir"),
  data.table(SEQN = done_ids("swan", cyc), cycle = cyc, tool = "swan")))))
avail <- avail[SEQN %in% an]
seqns <- sort(unique(avail$SEQN)); if (max_n > 0) seqns <- head(seqns, max_n)
cat("participants with at least one derived series:", length(seqns), " (GGIR", avail[tool == "ggir", .N], ", SWaN", avail[tool == "swan", .N], ")\n")

# ---- rebuilt series and PAXMIN rows of these participants -------------------------------------------------------------------------
ser <- as.data.table(readRDS(file.path(base, "01_recompute_SRI/intermediate/series_combined.rds")))[SEQN %in% seqns, .(SEQN, PAXDAYM, PAXPREDM, binary_sleep)]
ser[, PAXDAYM := as.integer(PAXDAYM)]
pax <- rbindlist(lapply(unique(avail[SEQN %in% seqns, cycle]), function(cyc) {
  x <- readRDS(paxrds[[cyc]]); setDT(x)
  x <- x[SEQN %in% seqns, .(SEQN, PAXDAYM = as.integer(PAXDAYM), PAXSSNMP, PAXTSM, PAXQFM, PAXAISMM, PAXPREDM = as.character(PAXPREDM))]
  gc(); x
}))
setorder(pax, SEQN, PAXSSNMP); setorder(ser, SEQN)
cat("rebuilt rows:", nrow(ser), " PAXMIN rows:", nrow(pax), "\n")

# ---- helpers ------------------------------------------------------------------------------------------------------------------------
sri_hand <- function(s, M = 1440L) {
  n <- length(s); if (n <= M) return(c(NA_real_, 0))
  a <- s[1:(n - M)]; b <- s[(M + 1):n]; ok <- !is.na(a) & !is.na(b); np <- sum(ok)
  if (np == 0) c(NA_real_, 0) else c(-100 + 200 * mean(a[ok] == b[ok]), np)
}
read_meta <- function(f) { x <- fread(f, header = FALSE, sep = "\t"); setNames(x$V2, x$V1) }
lab2bin <- function(x) fifelse(x == "sleep", 1L, fifelse(x == "wake", 0L, NA_integer_))
# GGIR: epochs of 5 s -> minute labels (majority; tie = wake; more than half of the epochs non-wear = missing)
ggir_minutes <- function(seqn, cyc) {
  d <- file.path(raw, "_ggir", cyc); mf <- file.path(d, "done", paste0(seqn, ".txt")); ef <- file.path(d, seqn, paste0(seqn, "_epochs.csv.gz"))
  if (!file.exists(mf) || !file.exists(ef)) return(NULL)
  m <- read_meta(mf); g <- fread(ef); codes <- fread(file.path(d, seqn, "behavioralcodes.csv"))
  g[, cn := codes$class_name[match(class_id, codes$class_id)]]
  t0 <- as.numeric(as.POSIXct(m[["raw_start"]], tz = "UTC")); g[, minute := as.integer(floor((timenum - t0) / 60))]
  g[, st := fifelse(invalidepoch == 1, "non-wear", fifelse(cn %in% c("spt_sleep", "day_sleep"), "sleep", "wake"))]
  mm <- g[, .(n = .N, sl = sum(st == "sleep"), nw = sum(st == "non-wear"), spt = mean(SleepPeriodTime)), by = minute]
  mm[, `:=`(GGIR_SIB = lab2bin(fifelse(nw > n / 2, "na", fifelse(sl > n / 2, "sleep", "wake"))),
            GGIR_SPT = lab2bin(fifelse(nw > n / 2, "na", fifelse(spt > 0.5, "sleep", "wake"))))]
  mm[, .(minute, GGIR_SIB, GGIR_SPT)]
}
# SWaN: two 30-s labels per minute; equal labels form the minute label, otherwise unknown (missing)
swan_minutes <- function(seqn, cyc) {
  f <- file.path(raw, "_swan", cyc, paste0(seqn, "_final.csv.gz")); if (!file.exists(f) || !file.exists(file.path(raw, "_swan", cyc, "done", paste0(seqn, ".txt")))) return(NULL)
  s <- fread(f); s[, START_TIME := as.POSIXct(START_TIME, tz = "UTC")]
  t0 <- min(s$START_TIME); s[, sec := as.numeric(difftime(START_TIME, t0, units = "secs"))]
  s[, `:=`(minute = as.integer(floor(sec / 60)), half = as.integer(round(sec %% 60) / 30))]
  s[, lab := c(Wear = "wake", Sleep = "sleep", Nonwear = "non-wear")[PREDICTION]]
  w <- dcast(s[, .(minute, half, lab)], minute ~ half, value.var = "lab"); setnames(w, c("minute", "w1", "w2"))
  w[, SWaN := lab2bin(fifelse(!is.na(w1) & !is.na(w2) & w1 == w2, w1, "unknown"))]
  w[, .(minute, SWaN)]
}
# fragmentation on the retained days: transitions per day on the observed sequence (missing minutes skipped), runs with missing as own state
frag_stats <- function(b, day) {
  d <- data.table(b = b, day = day)[!is.na(day)]
  tr <- d[, { s <- b[!is.na(b)]; .(tr = sum(diff(s) != 0), nv = sum(!is.na(b)), sl = sum(b == 1, na.rm = TRUE)) }, by = day][nv > 0]
  runs <- d[, { r <- rle(fifelse(is.na(b), -1L, b)); .(state = r$values, len = r$lengths) }, by = day]
  list(trans_day = mean(tr$tr), sleep_min_day = mean(tr$sl), runs = runs[state >= 0])
}

# ---- per participant ----------------------------------------------------------------------------------------------------------------
res <- list(); bouts <- list(); flags <- list()
serk <- split(ser, by = "SEQN", keep.by = FALSE); paxk <- split(pax, by = "SEQN", keep.by = FALSE)
for (k in seq_along(seqns)) {
  sq <- seqns[k]; cyc <- avail[SEQN == sq, cycle][1]
  r <- tryCatch({
    s <- serk[[as.character(sq)]]; p <- paxk[[as.character(sq)]]
    if (is.null(s) || is.null(p)) stop("no rebuilt series or PAXMIN rows")
    p <- p[PAXTSM == 60 & PAXQFM == 0 & PAXDAYM %in% unique(s$PAXDAYM)]
    if (nrow(p) != nrow(s)) stop(sprintf("row count differs: PAXMIN %d vs rebuilt %d", nrow(p), nrow(s)))
    if (!identical(p$PAXPREDM, as.character(s$PAXPREDM)) || !identical(p$PAXDAYM, s$PAXDAYM)) stop("labels or days differ between PAXMIN and the rebuilt series")
    p[, `:=`(m = as.integer(PAXSSNMP / 4800), pax = s$binary_sleep)]
    if (anyDuplicated(p$m)) stop("duplicated minute index")
    m0 <- min(p$m); L <- max(p$m) - m0 + 1L; ix <- p$m - m0 + 1L
    day <- rep(NA_integer_, L); day[ix] <- p$PAXDAYM
    pax_v <- rep(NA_integer_, L); pax_v[ix] <- as.integer(p$pax)
    ism_v <- rep(NA_real_, L); ism_v[ix] <- p$PAXAISMM
    ser_list <- list(); gm <- if (any(avail[SEQN == sq, tool] == "ggir")) ggir_minutes(sq, cyc) else NULL; sm <- if (any(avail[SEQN == sq, tool] == "swan")) swan_minutes(sq, cyc) else NULL
    put <- function(dt, col) { v <- rep(NA_integer_, L); j <- dt$minute - m0 + 1L; ok <- j >= 1L & j <= L; v[j[ok]] <- dt[[col]][ok]; v }
    if (!is.null(gm)) { ser_list$GGIR_SIB <- put(gm, "GGIR_SIB"); ser_list$GGIR_SPT <- put(gm, "GGIR_SPT") }
    if (!is.null(sm)) ser_list$SWaN <- put(sm, "SWaN")
    out <- list(); fr <- list()
    base_row <- data.table(SEQN = sq, cycle = cyc, n_days = uniqueN(p$PAXDAYM), n_ret = nrow(p), ism_share = mean(p$PAXAISMM) / 4800)
    sp <- sri_hand(pax_v); f0 <- frag_stats(pax_v, day)
    out$PAXPREDM <- cbind(base_row, deriv = "PAXPREDM", sri_A = sp[1], pairs_A = sp[2], valid_A = sum(!is.na(pax_v)), miss_in_mask = 0,
                          trans_day_A = f0$trans_day, sleep_min_day_A = f0$sleep_min_day, sri_B = sp[1], pairs_B = sp[2], miss_B = mean(is.na(pax_v[!is.na(day)])), incl_B = TRUE)
    fr$PAXPREDM <- f0$runs
    for (dn in names(ser_list)) {
      dv <- ser_list[[dn]]
      a <- fifelse(!is.na(pax_v) & !is.na(dv), dv, NA_integer_)                       # common mask
      b <- fifelse(!is.na(day), dv, NA_integer_)                                      # own missing minutes on the retained days
      sa <- sri_hand(a); sb <- sri_hand(b); fa <- frag_stats(a, day); missB <- mean(is.na(b[!is.na(day)]))
      ret <- !is.na(day) & !is.na(pax_v); ism <- ret & !is.na(ism_v) & ism_v >= 2400  # agreement among minutes that PAXPREDM labels wake or sleep
      cnt <- function(sel) c(N = sum(sel), agree = sum(sel & !is.na(dv) & dv == pax_v), dmiss = sum(sel & is.na(dv)),
                             n11 = sum(sel & dv == 1 & pax_v == 1, na.rm = TRUE), n10 = sum(sel & dv == 1 & pax_v == 0, na.rm = TRUE),
                             n01 = sum(sel & dv == 0 & pax_v == 1, na.rm = TRUE), n00 = sum(sel & dv == 0 & pax_v == 0, na.rm = TRUE))
      ca <- cnt(ret); ci <- cnt(ism); names(ca) <- paste0("all_", names(ca)); names(ci) <- paste0("ism_", names(ci))
      out[[dn]] <- cbind(base_row, deriv = dn, sri_A = sa[1], pairs_A = sa[2], valid_A = sum(!is.na(a)), miss_in_mask = ca[["all_dmiss"]] / ca[["all_N"]],
                         trans_day_A = fa$trans_day, sleep_min_day_A = fa$sleep_min_day, sri_B = sb[1], pairs_B = sb[2], miss_B = missB, incl_B = missB < 0.30,
                         as.data.table(as.list(c(ca, ci))))
      fr[[dn]] <- fa$runs
    }
    list(out = rbindlist(out, fill = TRUE), fr = rbindlist(lapply(names(fr), function(n) cbind(SEQN = sq, deriv = n, fr[[n]])))[, .(SEQN, deriv, state, len)])
  }, error = function(e) { flags[[length(flags) + 1]] <<- data.table(SEQN = sq, reason = conditionMessage(e)); NULL })
  if (!is.null(r)) { res[[length(res) + 1]] <- r$out; bouts[[length(bouts) + 1]] <- r$fr }
  if (k %% 50 == 0) cat(sprintf("  %d/%d (%.1f min)\n", k, length(seqns), as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
}
res <- rbindlist(res, fill = TRUE); bouts <- rbindlist(bouts); flags <- if (length(flags)) rbindlist(flags) else data.table(SEQN = integer(), reason = character())
saveRDS(res, file.path(outdir, "raw_derivation_participants.rds")); fwrite(res, file.path(outdir, "raw_derivation_participants.csv"))
saveRDS(bouts, file.path(outdir, "raw_derivation_bouts.rds")); fwrite(flags, file.path(outdir, "raw_derivation_flags.csv"))
cat("done:", uniqueN(res$SEQN), "participants,", nrow(res), "rows; flagged:", nrow(flags), "\n")
print(res[, .(n = .N, sri_A = round(mean(sri_A, na.rm = TRUE), 1), miss_in_mask = round(mean(miss_in_mask, na.rm = TRUE), 3), trans_day = round(mean(trans_day_A, na.rm = TRUE), 1)), by = deriv])
# computing environment
sd <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (!is.na(sd)) try(system2("Rscript", c(file.path(dirname(sd), "system_info.R"), file.path(outdir, "system_info_raw_derivation.txt")), stdout = FALSE), silent = TRUE)
