suppressMessages(library(GGIR))
# GGIR's own unit-test scenarios (tests/testthat/test_CalcSleepRegularityIndex.R, GitHub main), fed through OUR wrapper call style
tz <- "Europe/Amsterdam"
time <- seq(from = as.POSIXct("2020-01-01 00:01:10", tz = tz), to = as.POSIXct("2020-01-06 00:01:10", tz = tz), by = 5)
NID <- 12*1440; sleep <- rep(0, length(time)); R1 <- round(NID*0.4); R2 <- round(NID*0.5); R3 <- round(NID*0.8); R4 <- round(NID*0.9)
sleep[R1:R3] <- 1; sleep[(NID+R2):(NID+R3)] <- 1; sleep[(NID*2+R2):(NID*2+R4)] <- 1; sleep[(NID*3+R1):(NID*3+R4)] <- 1; sleep[(NID*4+R1):(NID*4+R3)] <- 1; sleep[(NID*5+R2):(NID*5+R3)] <- 1
sleep <- sleep[1:length(time)]
our_call <- function(ts, s, epochsize, tz) { df <- data.frame(time = ts, sleepstate = as.numeric(s), invalid = as.integer(is.na(s)))
  GGIR::CalcSleepRegularityIndex(data = df, epochsize = epochsize, desiredtz = tz) }
S1 <- our_call(time, sleep, 5, tz); cat("Scenario 1 (expected day-pairs 2 and 3 = 80):\n"); print(S1[, c("day","SleepRegularityIndex","frac_valid")])
S2 <- our_call(time, rep(0, length(time)), 5, tz); cat("Scenario 2 (expected 100):", S2$SleepRegularityIndex[3:4], "\n")
S3 <- our_call(time, 1:length(time), 5, tz); cat("Scenario 3 (expected -100):", S3$SleepRegularityIndex[2:3], "\n")
# Hand Phillips on the same series (position-based, 5-s epochs, M = 17280) for day-pairs 2,3 equivalents
hand_pair <- function(s, M) { i <- 1:(length(s)-M); a <- s[i]; b <- s[i+M]; v <- !is.na(a) & !is.na(b); -100 + 200*mean(a[v]==b[v]) }
cat("Hand Phillips over whole series (5-s):", round(hand_pair(sleep, 17280), 2), " | mean of GGIR day-pairs:", round(mean(S1$SleepRegularityIndex, na.rm=TRUE), 2), "\n")
# Minute-level version of the same test (our NHANES usage: epochsize = 60)
tm <- seq(from = as.POSIXct("2020-01-01 00:00:00", tz = "UTC"), by = 60, length.out = 6*1440)
sm <- rep(0, length(tm)); D <- 1440; r1 <- round(D*0.4); r2 <- round(D*0.5); r3 <- round(D*0.8); r4 <- round(D*0.9)
sm[r1:r3] <- 1; sm[(D+r2):(D+r3)] <- 1; sm[(2*D+r2):(2*D+r4)] <- 1; sm[(3*D+r1):(3*D+r4)] <- 1; sm[(4*D+r1):(4*D+r3)] <- 1; sm[(5*D+r2):(5*D+r3)] <- 1
G <- our_call(tm, sm, 60, "UTC"); cat("Minute-level analogue, GGIR day-pairs:", round(G$SleepRegularityIndex, 2), "\n")
cat("Minute-level analogue, hand per day-pair:", sapply(1:5, function(d) round(hand_pair(sm[((d-1)*D+1):((d+1)*D)], D), 2)), "\n")
# Missing calendar day: drop day 3 entirely -> what does GGIR pair?
keep <- !(seq_along(tm) %in% ((2*D+1):(3*D)))
Gm <- our_call(tm[keep], sm[keep], 60, "UTC"); cat("With calendar day 3 removed, GGIR day-pairs:\n"); print(Gm[, c("day","SleepRegularityIndex","frac_valid","date")])
