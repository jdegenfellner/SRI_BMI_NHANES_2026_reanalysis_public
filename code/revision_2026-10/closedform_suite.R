# Closed-form checks of the three calculators (hand / GGIR / sleepreg) on deterministic minute series with known SRI,
# each presented (A) as consecutive full days and (B) with one middle day marked NA. Days are midnight-aligned; the
# recording starts at 12:00 on day 0 (placement) with NA until midnight, like the NHANES NA timeline.
suppressMessages({library(data.table); library(GGIR)})
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())   # working tree root (contains 00_original_paper_materials/ and 01_recompute_SRI/)
source(file.path(base, "01_recompute_SRI/sleepreg_src/sleepreg_allFunctions.R")); load(file.path(base, "01_recompute_SRI/sleepreg_src/sysdata.rda"))
D <- 1440L
day <- function(on_min, off_min, nap = NULL) { s <- integer(D); s[((on_min:(D-1)) %% D) + 1] <- 1L; s[1:off_min] <- 1L; if (!is.null(nap)) s[nap[1]:nap[2] + 1] <- 1L; s }   # sleep from on_min (evening) through midnight to off_min (morning)
series <- function(days) c(rep(NA_integer_, 720), unlist(days))   # placement at noon of day 0, NA until midnight
ts_of <- function(s) 12*3600 + seq(0, by = 60, length.out = length(s))
hand <- function(s, M = D) { n <- length(s); i <- 1:(n-M); a <- s[i]; b <- s[i+M]; ok <- !is.na(a) & !is.na(b); -100 + 200*mean(a[ok] == b[ok]) }
ggir <- function(ts, s) { df <- data.frame(time = as.POSIXct(ts, origin="1970-01-01", tz="UTC"), sleepstate = as.numeric(s), invalid = as.integer(is.na(s))); o <- GGIR::CalcSleepRegularityIndex(data = df, epochsize = 60, desiredtz = "UTC"); mean(o$SleepRegularityIndex[o$frac_valid >= 0.66]) }
sleepreg <- function(ts, s, tag) { bd <- file.path(getwd(), paste0("cf_", tag)); unlink(bd, recursive=TRUE); dir.create(bd); fwrite(data.table(b = s, t = ts), file.path(bd, "sw_1.csv"))
  invisible(capture.output({ SWS_from_binarySW(binarySWdir = bd, tsCol = 2, binaryCol = 1); od <- file.path(bd, "out"); dir.create(od); SRI_from_binary(binarydir = file.path(bd, "SWS_output"), alloutdir = od, col.trans = 1, col.timestamp = 2, wr.raster = FALSE, overwr = TRUE, minSRIdays = 3) }))
  as.numeric(fread(file.path(od, "SRI.csv"))$SRI[1]) }
run <- function(label, days, expected) {
  sA <- series(days); tA <- ts_of(sA)
  sB <- sA; sB[(720 + 3*D + 1):(720 + 4*D)] <- NA   # day 4 invalid
  rA <- c(hand(sA), ggir(tA, sA), sleepreg(tA, sA, "A")); rB <- c(hand(sB), ggir(tA, sB), sleepreg(tA, sB, "B"))
  cat(sprintf("%-46s expected %7.3f | A: hand %7.3f GGIR %7.3f sleepreg %7.3f | B (day 4 NA): hand %7.3f GGIR %7.3f sleepreg %7.3f\n", label, expected, rA[1], rA[2], rA[3], rB[1], rB[2], rB[3])) }
# 1. identical days (sleep 23:00-07:00): SRI = 100
run("1 identical days", replicate(7, day(23*60, 7*60), simplify = FALSE), 100)
# 2. alternating 30-min shift (sleep onset/offset shift by +30 on odd days): each pair disagrees for 2*30 min -> 100 - 200*60/1440 = 91.667
run("2 alternating +-30 min shift of the whole block", lapply(1:7, function(d) day(23*60 + 30*(d %% 2), 7*60 + 30*(d %% 2))), 100 - 200*60/1440)
# 3. constant drift +15 min/day (onset and offset): consecutive days disagree 30 min -> 100 - 200*30/1440 = 95.833
run("3 drift +15 min/day (onset 22:00 ->23:30)", lapply(1:7, function(d) day(22*60 + 15*(d-1), 6*60 + 15*(d-1))), 100 - 200*30/1440)
# 4. nap 14:00-14:30 on days 1,3,5,7 only: pairs (1,2),(2,3),...: every pair disagrees 30 min -> 100 - 200*30/1440 = 95.833
run("4 30-min nap on alternate days", lapply(1:7, function(d) day(23*60, 7*60, if (d %% 2 == 1) c(14*60, 14*60+29))), 100 - 200*30/1440)
# 5. perfectly inverted alternate days (sleep 23-07 vs sleep 11-19): each pair disagrees in 2*8 h = 960 min -> 100 - 200*960/1440 = -33.33
run("5 night sleep vs day sleep on alternate days", lapply(1:7, function(d) if (d %% 2 == 1) day(23*60, 7*60) else { s <- integer(D); s[(11*60):(19*60-1) + 1] <- 1L; s }), 100 - 200*960/1440)
# 6. random pattern each day (p = 0.5): expected 0 (sampling noise ~ +-1)
set.seed(7); run("6 random 50% sleep (expected ~0)", replicate(7, as.integer(rbinom(D, 1, 0.5)), simplify = FALSE), 0)
