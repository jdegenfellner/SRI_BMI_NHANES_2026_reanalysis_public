suppressMessages({library(data.table); library(survey)})
options(width = 170)
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())   # working tree root (contains 00_original_paper_materials/ and 01_recompute_SRI/)
origdata <- file.path(base, "00_original_paper_materials", "data")
d  <- as.data.table(readRDS(file.path(origdata, "df_an1_imputed_MAIN.RDS")))
an <- as.data.table(readRDS(file.path(origdata, "Combined_data_sets_2_waves_analytic_data_set.RDS")))
d  <- merge(d, an[, .(SEQN, WTMEC2YR, SDMVPSU, SDMVSTRA, SRI_orig = SRI_no_imputation)], by = "SEQN", all.x = TRUE)[, WTMEC4YR := WTMEC2YR/2]
an_seqn <- unique(d$SEQN)

ser <- as.data.table(readRDS(file.path(base, "01_recompute_SRI/intermediate/series_combined.rds")))[SEQN %in% an_seqn]
setorder(ser, SEQN, UNIX_TIMESTAMP)
ser[, dayi := as.integer(PAXDAYM)]
ser[, day_start := UNIX_TIMESTAMP[1], by = .(SEQN, dayi)]
ser[, tod_min := as.integer(round((UNIX_TIMESTAMP - day_start)/60))]
# contiguous blocks of consecutive calendar days within a participant
ser[, blk := cumsum(c(1L, as.integer(diff(dayi) > 1L))), by = SEQN]   # new block only when a calendar day is missing
ser[, noon_day := fifelse(tod_min >= 720L, dayi, dayi - 1L)]
cat("participants:", uniqueN(ser$SEQN), " rows:", nrow(ser), " blocks:", uniqueN(ser[, .(SEQN, blk)]), "\n")

# ---------- pass-2 operators: NA-transparent bridging, majority filter, GGIR exact rule ----------
roll_frac_sleep <- function(x, w = 60L) { n <- length(x); isn <- !is.na(x); xs <- ifelse(isn, x, 0)
  cs <- c(0, cumsum(xs)); cn <- c(0, cumsum(isn)); h <- w %/% 2L
  lo <- pmax(1L, seq_len(n) - h); hi <- pmin(n, seq_len(n) + h); s <- cs[hi + 1L] - cs[lo]; m <- cn[hi + 1L] - cn[lo]; ifelse(m > 0, s/m, NA_real_) }
v_thresh <- function(x, w = 60L, frac) { f <- roll_frac_sleep(x, w); y <- ifelse(is.na(x) | is.na(f), NA, ifelse(f >= frac, 1, 0)); y }   # GGIR-style: >=frac -> sleep, else wake
na_fill <- function(x, L = 15L) {   # NA runs <= L min with identical neighbours take that state
  r <- rle(ifelse(is.na(x), -1, x)); v <- r$values; l <- r$lengths; k <- length(v)
  if (k >= 3) for (i in 2:(k-1)) if (v[i] == -1 && l[i] <= L && v[i-1] == v[i+1] && v[i-1] != -1) v[i] <- v[i-1]
  y <- inverse.rle(list(values = v, lengths = l)); y[y == -1] <- NA; y }
v_close <- function(x, W = 30L) { r <- rle(ifelse(is.na(x), -1, x)); v <- r$values; l <- r$lengths; k <- length(v)
  if (k >= 3) for (i in 2:(k-1)) if (v[i] == 0 && l[i] < W && v[i-1] == 1 && v[i+1] == 1) v[i] <- 1
  y <- inverse.rle(list(values = v, lengths = l)); y[y == -1] <- NA; y }
v_open <- function(x, S = 30L) { r <- rle(ifelse(is.na(x), -1, x)); v <- r$values; l <- r$lengths; k <- length(v)
  if (k >= 3) for (i in 2:(k-1)) if (v[i] == 1 && l[i] < S && v[i-1] == 0 && v[i+1] == 0) v[i] <- 0
  y <- inverse.rle(list(values = v, lengths = l)); y[y == -1] <- NA; y }
v_bridgeNA <- function(x, W = 30L, S = 30L) { y <- na_fill(x); y <- v_close(y, W); y <- na_fill(y); y <- v_close(y, W); y <- na_fill(y); y <- v_open(y, S); y <- na_fill(y); v_open(y, S) }
v_bridge2 <- function(x, W = 30L, S = 30L) {   # bridging on the observed (non-NA) sequence, so unknown minutes do not block flanking
  idx <- which(!is.na(x)); if (length(idx) < 3) return(x)
  z <- x[idx]; r <- rle(z); v <- r$values; l <- r$lengths; k <- length(v)
  if (k >= 3) for (i in 2:(k-1)) if (v[i] == 0 && l[i] < W && v[i-1] == 1 && v[i+1] == 1) v[i] <- 1
  z <- inverse.rle(list(values = v, lengths = l)); r <- rle(z); v <- r$values; l <- r$lengths; k <- length(v)
  if (k >= 3) for (i in 2:(k-1)) if (v[i] == 1 && l[i] < S && v[i-1] == 0 && v[i+1] == 0) v[i] <- 0
  y <- x; y[idx] <- inverse.rle(list(values = v, lengths = l)); na_fill(y, 15L) }
v_spt_only <- function(x, noon_day) { r <- rle(ifelse(is.na(x), -1, x)); v <- r$values; l <- r$lengths
  ends <- cumsum(l); starts <- ends - l + 1L; nd <- noon_day[starts]; keep <- rep(FALSE, length(v)); is1 <- which(v == 1)
  if (length(is1)) { dt <- data.table(i = is1, nd = nd[is1], l = l[is1]); keep[dt[, i[which.max(l)], by = nd]$V1] <- TRUE }
  v[v == 1 & !keep] <- 0; y <- inverse.rle(list(values = v, lengths = l)); y[y == -1] <- NA; y }
# ---------- hand SRI (replicates the published calculation: weekday-consecutive pairs, 1-min epochs) ----------
hand_sri <- function(s, dow, M = 1440L) {
  n <- length(s); if (n <= M) return(NA_real_)
  i <- 1:(n - M); cur <- dow[i]; nxt <- dow[i + M]
  consec <- (nxt - cur == 1) | (cur == 7 & nxt == 1)
  a <- s[i]; b <- s[i + M]; valid <- consec & !is.na(a) & !is.na(b)
  Nv <- sum(valid); if (Nv == 0) return(NA_real_)
  -100 + 200/Nv * sum(a[valid] == b[valid])
}
frag_stats <- function(x) {                         # per contiguous vector: transitions/day, longest sleep run, sleep min/day
  r <- rle(ifelse(is.na(x), -1, x)); v <- r$values; l <- r$lengths
  s <- v[v != -1]
  c(trans = sum(diff(s) != 0), longest = if (any(v == 1)) max(l[v == 1]) else 0, sleep = sum(l[v == 1]), ndays = length(x)/1440)
}

# ---------- build the variants per participant (per contiguous block for smoothing) ----------
ser[, raw := binary_sleep]
ser[, maj1h50 := v_thresh(raw, 60L, 0.5), by = .(SEQN, blk)]
ser[, ggir1h80 := v_thresh(raw, 60L, 0.8), by = .(SEQN, blk)]
ser[, bridge2 := v_bridge2(raw, 30L, 30L), by = .(SEQN, blk)]
ser[, spt2 := v_spt_only(bridge2, noon_day), by = .(SEQN, blk)]
variants <- c("raw", "maj1h50", "ggir1h80", "bridge2", "spt2")
ser[, ts_fix := (dayi - 1) * 86400 + tod_min * 60]
saveRDS(ser[, c("SEQN", "dayi", "PAXDAYWM", "ts_fix", variants), with = FALSE], "series_variants.rds")
cat("series_variants.rds saved\n")
cat("variants built\n")

res <- ser[, {
  out <- list()
  for (v in variants) {
    x <- get(v)
    out[[paste0("SRI_", v)]] <- hand_sri(x, PAXDAYWM)
    fs <- frag_stats(x)
    out[[paste0("trans_", v)]] <- fs["trans"]/fs["ndays"]
    out[[paste0("longest_", v)]] <- fs["longest"]
    out[[paste0("sleepmin_", v)]] <- fs["sleep"]/fs["ndays"]
  }
  out
}, by = SEQN]
res <- merge(res, d[, .(SEQN, SRI_orig)], by = "SEQN")
saveRDS(res, "variants4_sri.rds"); fwrite(res, "variants4_sri.csv")

cat("\n==== validation: hand SRI on raw vs published SRI_orig ====\n")
cat("r =", round(cor(res$SRI_raw, res$SRI_orig, use = "pair"), 5), " max|diff| =", round(max(abs(res$SRI_raw - res$SRI_orig), na.rm = TRUE), 4), "\n")
cat("\n==== SRI distribution by input variant (n =", nrow(res), ") ====\n")
for (v in variants) {
  x <- res[[paste0("SRI_", v)]]
  cat(sprintf("%-9s median %5.1f  mean %5.1f  sd %4.1f  Q1 %5.1f  Q3 %5.1f  min %6.1f  | r with raw %.3f | transitions/day %4.1f | longest sleep run (median min) %4.0f | sleep min/day %4.0f\n",
      v, median(x, na.rm=TRUE), mean(x, na.rm=TRUE), sd(x, na.rm=TRUE), quantile(x,.25,na.rm=TRUE), quantile(x,.75,na.rm=TRUE), min(x,na.rm=TRUE),
      cor(x, res$SRI_raw, use="pair"), mean(res[[paste0("trans_",v)]]), median(res[[paste0("longest_",v)]]), mean(res[[paste0("sleepmin_",v)]])))
}

# ---------- Table-3-style models (survey-weighted, published recipe) ----------
set.seed(3455); SRIsim <- replicate(100000, -100 + 200*sum(rbinom(1440,1,.5))/1440)
m_s <- mean(SRIsim); s_s <- sd(SRIsim)
winz <- function(x){ i <- which(x < -6.8); x[i] <- rnorm(length(i), m_s, s_s); x }
cov_full <- c("RIDAGEYR","RIAGENDR_factor","RIDRETH3_factor","DMDEDUC2_factor","INDHHIN2_factor",
              "OCD150_factor","alcohol_consumption","smoking_status","LBXVIDMS","DR1TKCAL_winsorized",
              "Depression_score_PHQ_9","DMDMARTL_factor","activity_level")
cov_full_nosex <- setdiff(cov_full, "RIAGENDR_factor")
beta_SRI <- function(des, rhs) coef(svyglm(as.formula(paste0("log(BMXBMI) ~ SRI + ", rhs)), design = des))
fit_set <- function(dd) {
  des <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC4YR, nest=TRUE, data=dd)
  b1 <- beta_SRI(des, "RIDAGEYR + RIAGENDR_factor")["SRI"]
  b2 <- beta_SRI(des, paste(cov_full, collapse=" + "))["SRI"]
  c3 <- beta_SRI(des, paste0("RIAGENDR_factor*SRI - SRI + SRI + ", paste(cov_full_nosex, collapse=" + ")))
  ix <- grep("RIAGENDR_factor.*:SRI|SRI:RIAGENDR_factor", names(c3), value=TRUE)
  bref <- c3["SRI"]; bother <- bref + (if (length(ix)) c3[ix[1]] else 0)
  c(b1 = unname(b1), b2 = unname(b2), bref = unname(bref), bother = unname(bother))
}
mf_table <- function(dd, svals, nb = 200, seed = 11) {
  dd <- copy(dd); dd[, SRI := winz(svals)]
  dd <- dd[is.finite(SRI) & !is.na(BMXBMI) & !is.na(WTMEC4YR)]
  med <- dd[, median(SRI), by=.(q=cut(SRI, quantile(SRI,0:5/5), include.lowest=TRUE, labels=1:5))][order(q), V1]
  dSRI <- med[5] - med[1]
  pt <- exp(fit_set(dd) * dSRI)
  set.seed(seed); n <- nrow(dd)
  B <- t(replicate(nb, { bs <- dd[sample(n, replace=TRUE)]; tryCatch(exp(fit_set(bs) * dSRI), error=function(e) rep(NA,4)) }))
  ci <- apply(B, 2, quantile, c(.025,.975), na.rm=TRUE)
  reflev <- levels(dd$RIAGENDR_factor)[1]; otherlev <- levels(dd$RIAGENDR_factor)[2]
  data.table(n = n, dSRI = round(dSRI,1), medQ1 = round(med[1],1), medQ5 = round(med[5],1),
             M1 = sprintf("%.3f (%.3f, %.3f)", pt[1], ci[1,1], ci[2,1]),
             M2 = sprintf("%.3f (%.3f, %.3f)", pt[2], ci[1,2], ci[2,2]),
             M3_ref = sprintf("%s %.3f (%.3f, %.3f)", reflev, pt[3], ci[1,3], ci[2,3]),
             M3_other = sprintf("%s %.3f (%.3f, %.3f)", otherlev, pt[4], ci[1,4], ci[2,4]))
}
dm <- merge(d, res, by = "SEQN")
rows <- list()
rows[["published"]] <- cbind(variant = "published (orig, winsorised)", mf_table(dm, dm$SRI_no_imputation_winsorized_rand)); cat("models: published done\n")
for (v in variants) { rows[[v]] <- cbind(variant = v, mf_table(dm, dm[[paste0("SRI_", v)]])); cat("models:", v, "done\n") }
tab <- rbindlist(rows)
fwrite(tab, "variants4_models.csv"); saveRDS(tab, "variants4_models.rds")
cat("\n==== MF Q1->Q5 on BMI by sleep/wake INPUT variant (hand/Phillips SRI, published Table-3 recipe, 200 bootstrap) ====\n")
print(tab)
cat("\nDONE\n")
