# Gap presentation test, midnight-aligned: placement at noon of day 0 (NA until midnight), 7 full days, end at midnight.
suppressMessages({library(data.table); library(GGIR)})
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())   # working tree root (contains 00_original_paper_materials/ and 01_recompute_SRI/)
source(file.path(base, "01_recompute_SRI/sleepreg_src/sleepreg_allFunctions.R")); load(file.path(base, "01_recompute_SRI/sleepreg_src/sysdata.rda"))
D <- 1440L; shifts <- c(0, 30, -45, 15, 60, -20, 10); naps <- c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE); s7 <- integer(7*D)
for (d in 1:7) { on <- (d-1)*D + 23*60+30 + shifts[d]; off <- d*D + 7*60 + shifts[d]; idx <- on:(off-1); idx <- idx[idx >= 0 & idx < 7*D]; s7[idx + 1] <- 1L; if (naps[d]) { ni <- (d-1)*D + 14*60 + 0:19; s7[ni + 1] <- 1L } }
sA <- s7; tA <- 86400 + seq(0, by = 60, length.out = length(sA))   # ends absent: series starts 00:00 of day 1
hand <- function(s, M = D) { n <- length(s); i <- 1:(n-M); a <- s[i]; b <- s[i+M]; ok <- !is.na(a) & !is.na(b); -100 + 200*mean(a[ok] == b[ok]) }
ggir <- function(ts, s) { df <- data.frame(time = as.POSIXct(ts, origin="1970-01-01", tz="UTC"), sleepstate = as.numeric(s), invalid = as.integer(is.na(s))); GGIR::CalcSleepRegularityIndex(data = df, epochsize = 60, desiredtz = "UTC") }
sleepreg <- function(ts, s, tag, md) { ts <- ts + 30; bd <- file.path(getwd(), paste0("g2_", tag)); unlink(bd, recursive=TRUE); dir.create(bd); fwrite(data.table(b = s, t = ts), file.path(bd, "sw_1.csv"))
  invisible(capture.output({ SWS_from_binarySW(binarySWdir = bd, tsCol = 2, binaryCol = 1); od <- file.path(bd, "out"); dir.create(od); SRI_from_binary(binarydir = file.path(bd, "SWS_output"), alloutdir = od, col.trans = 1, col.timestamp = 2, wr.raster = FALSE, overwr = TRUE, minSRIdays = md) }))
  r <- fread(file.path(od, "SRI.csv")); sprintf("%.2f (days %.0f)", as.numeric(r$SRI[1]), as.numeric(r$SRI_days[1])) }
show <- function(label, ts, s) { g <- ggir(ts, s); ok <- g$frac_valid >= 0.66
  cat(sprintf("%-32s hand %6.2f | GGIR filtered %6.2f  unfiltered %6.2f | sleepreg min3 %s  min2 %s\n   GGIR pairs: %s\n", label, hand(s), mean(g$SleepRegularityIndex[ok]), mean(g$SleepRegularityIndex), sleepreg(ts, s, paste0(gsub("[^a-z]", "", tolower(label)), 3), 3), sleepreg(ts, s, paste0(gsub("[^a-z]", "", tolower(label)), 2), 2), paste(sprintf("%s:%.1f(%.2f)", g$date, g$SleepRegularityIndex, g$frac_valid), collapse = " "))) }
show("complete 7 days", tA, sA)
hand_wd <- function(s, dow, M = D) { n <- length(s); i <- 1:(n-M); consec <- (dow[i+M]-dow[i]==1)|(dow[i]==7&dow[i+M]==1); a <- s[i]; b <- s[i+M]; v <- consec & !is.na(a) & !is.na(b); -100 + 200*mean(a[v]==b[v]) }
sB <- sA; sB[(3*D + 1):(4*D)] <- NA; show("day 4 as missing minutes", tA, sB)
rm4 <- !(seq_along(sA) %in% ((3*D + 1):(4*D))); show("day 4 deleted", tA[rm4], sA[rm4])
dow <- rep(c(3,4,5,6,7,1,2), each = D); cat("published weekday-rule formula on the deleted-day series:", round(hand_wd(sA[rm4], dow[rm4]), 2), "\n")
