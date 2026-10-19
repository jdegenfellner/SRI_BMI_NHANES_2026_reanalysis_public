suppressMessages(library(data.table)); D <- 1440L
day <- function(on_min, off_min, nap = NULL) { s <- integer(D); s[((on_min:(D-1)) %% D) + 1] <- 1L; s[1:off_min] <- 1L; if (!is.null(nap)) s[nap[1]:nap[2] + 1] <- 1L; s }
cases <- list(
  identical = replicate(7, day(23*60, 7*60), simplify = FALSE),
  shift30 = lapply(1:7, function(d) day(23*60 + 30*(d %% 2), 7*60 + 30*(d %% 2))),
  drift15 = lapply(1:7, function(d) day(22*60 + 15*(d-1), 6*60 + 15*(d-1))),
  nap = lapply(1:7, function(d) day(23*60, 7*60, if (d %% 2 == 1) c(14*60, 14*60+29))),
  inverted = lapply(1:7, function(d) if (d %% 2 == 1) day(23*60, 7*60) else { s <- integer(D); s[(11*60):(19*60-1) + 1] <- 1L; s }))
set.seed(7); cases$random <- replicate(7, as.integer(rbinom(D, 1, 0.5)), simplify = FALSE)
shifts <- c(0, 30, -45, 15, 60, -20, 10); naps <- c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE); s7 <- integer(7*D)
for (d in 1:7) { on <- (d-1)*D + 23*60+30 + shifts[d]; off <- d*D + 7*60 + shifts[d]; idx <- on:(off-1); idx <- idx[idx >= 0 & idx < 7*D]; s7[idx + 1] <- 1L; if (naps[d]) { ni <- (d-1)*D + 14*60 + 0:19; s7[ni + 1] <- 1L } }
out <- rbindlist(c(lapply(names(cases), function(n) data.table(case = n, minute = seq_len(7*D) - 1L, s = unlist(cases[[n]]))), list(data.table(case = "realistic", minute = seq_len(7*D) - 1L, s = s7))))
fwrite(out, "synthetic_cases.csv"); cat("exported", uniqueN(out$case), "cases\n")
hand <- function(s, M = D) { n <- length(s); i <- 1:(n-M); -100 + 200*mean(s[i] == s[i+M]) }
print(out[, .(hand_SRI = round(hand(s), 3)), by = case])
