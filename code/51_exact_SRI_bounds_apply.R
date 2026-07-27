# =====================================================================
# Multiverse paper — apply the VERIFIED exact SRI identification bounds
# to real NHANES participants. For each person compute [SRI_min, SRI_max]
# (sharp identification region under missingness) and compare to the
# published point SRI and the heuristic max/min imputations.
# Run AFTER part 2 finishes (shares the series_*.rds).
# =====================================================================
suppressMessages({library(data.table)})
options(width = 150)
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())
inter <- file.path(base, "01_recompute_SRI", "intermediate")

# verified vectorised multi-column DP (rows = consecutive days)
matrix_opt <- function(X, adj, want) {
  D <- nrow(X)
  INF  <- if (want=="max") -Inf else Inf
  pick <- if (want=="max") pmax else pmin
  init <- function(r){ f0<-rep(0,length(r)); f1<-rep(0,length(r))
                       f0[which(r==1)]<-INF; f1[which(r==0)]<-INF; list(f0=f0,f1=f1) }
  st <- init(X[1,]); total <- rep(0, ncol(X))
  for (d in 2:D) {
    if (!adj[d-1]) { total <- total + pick(st$f0, st$f1); st <- init(X[d,]); next }
    nf0 <- pick(st$f0+1, st$f1); nf1 <- pick(st$f1+1, st$f0)
    r <- X[d,]; nf0[which(r==1)] <- INF; nf1[which(r==0)] <- INF
    st <- list(f0=nf0, f1=nf1)
  }
  sum(total + pick(st$f0, st$f1))
}

# per-participant exact bounds
person_bounds <- function(s, paxdaym, paxdaywm) {
  D <- length(s)/1440
  X  <- t(matrix(s, nrow = 1440, ncol = D))          # D x 1440 (rows = days)
  dw <- as.numeric(paxdaywm[seq(1, length(paxdaywm), by = 1440)])  # weekday per day
  adj <- (diff(dw) == 1) | (head(dw,-1)==7 & tail(dw,-1)==1)        # consecutive-day adjacency
  P <- 1440 * sum(adj)                                # total comparable pairs (denominator)
  if (P == 0) return(c(NA, NA, NA))
  maxA <- matrix_opt(X, adj, "max"); minA <- matrix_opt(X, adj, "min")
  c(SRI_min = -100 + 200*minA/P, SRI_max = -100 + 200*maxA/P, pairs = P)
}

# ---- load series (analytic sample) ----
an <- as.data.table(readRDS(file.path(base,"00_original_paper_materials/data/df_an1_imputed_MAIN.RDS")))
anseq <- unique(an$SEQN)
s1 <- as.data.table(readRDS(file.path(inter,"series_1112.rds")))[SEQN %in% anseq]
s2 <- as.data.table(readRDS(file.path(inter,"series_1314.rds")))[SEQN %in% anseq]
ser <- rbind(s1,s2); rm(s1,s2); gc()
setorder(ser, SEQN, UNIX_TIMESTAMP); setkey(ser, SEQN)
seqns <- sort(unique(ser$SEQN))
cat("participants:", length(seqns), "\n")

res <- data.table(SEQN=seqns, SRI_min=NA_real_, SRI_max=NA_real_, pairs=NA_real_)
t0 <- Sys.time()
for (k in seq_along(seqns)) {
  d <- ser[.(seqns[k])]
  b <- person_bounds(d$binary_sleep, d$PAXDAYM, d$PAXDAYWM)
  res$SRI_min[k]<-b[1]; res$SRI_max[k]<-b[2]; res$pairs[k]<-b[3]
  if (k %% 1000 == 0) cat(sprintf("  %d/%d (%.1f min)\n", k, length(seqns),
        as.numeric(difftime(Sys.time(),t0,units="mins"))))
}
saveRDS(res, file.path(inter,"sri_exact_bounds.rds"))

# ---- compare to published point + heuristic max/min ----
ad <- as.data.table(readRDS(file.path(base,
  "00_original_paper_materials/data/Combined_data_sets_2_waves_analytic_data_set.RDS")))
m <- merge(res, ad[, .(SEQN, SRI_pt=SRI_no_imputation, SRI_heurmax=SRI_maximize, SRI_heurmin=SRI_minimize)],
           by="SEQN", all.x=TRUE)
m[, width := SRI_max - SRI_min]
cat("\n=== identification interval [SRI_min, SRI_max] ===\n")
cat(sprintf("median width: %.2f   IQR [%.2f, %.2f]   max width: %.2f\n",
    median(m$width,na.rm=T), quantile(m$width,.25,na.rm=T), quantile(m$width,.75,na.rm=T), max(m$width,na.rm=T)))
cat(sprintf("point SRI inside [min,max]: %.1f%%\n", 100*mean(m$SRI_pt>=m$SRI_min-1e-6 & m$SRI_pt<=m$SRI_max+1e-6, na.rm=T)))
cat("\n=== exactness check vs heuristic (must hold for ALL) ===\n")
cat("exact_max >= heuristic_max:", all(m$SRI_max >= m$SRI_heurmax-1e-6, na.rm=TRUE),
    " | exact_min <= heuristic_min:", all(m$SRI_min <= m$SRI_heurmin+1e-6, na.rm=TRUE), "\n")
cat(sprintf("mean gain exact_max - heuristic_max: %.3f ; exact_min - heuristic_min: %.3f\n",
    mean(m$SRI_max - m$SRI_heurmax, na.rm=T), mean(m$SRI_min - m$SRI_heurmin, na.rm=T)))
cat(sprintf("intervals crossing random-sleep 0: %.1f%% ; crossing winsor cut -6.8: %.1f%%\n",
    100*mean(m$SRI_min<0 & m$SRI_max>0, na.rm=T), 100*mean(m$SRI_min< -6.8 & m$SRI_max> -6.8, na.rm=T)))
fwrite(m, file.path(inter,"sri_exact_bounds_compare.csv"))
cat("\nDONE exact bounds application\n")
