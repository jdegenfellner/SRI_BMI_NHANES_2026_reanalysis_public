# =====================================================================
# Multiverse paper — METHOD CORE: exact SRI identification bounds.
# The SRI couples epochs only within the same time-of-day column across
# consecutive days => problem decomposes into independent binary chains.
# Each chain: maximize/minimize # equal consecutive pairs with some fixed
# entries, solved by a 2-state DP. THIS SCRIPT VERIFIES the DP against
# brute-force enumeration on small chains, then against a vectorised
# multi-column implementation. (No big data loaded.)
# =====================================================================
options(width = 140)
set.seed(1)

# ---- (a) single-chain DP: opt # of equal adjacent pairs ----
# x: vector of 0/1/NA (NA = free to impute). want="max" or "min".
chain_opt <- function(x, want = "max") {
  L <- length(x); if (L <= 1) return(0L)
  INF <- if (want == "max") -Inf else Inf
  cmb <- if (want == "max") max else min
  # f[1]=best agreements for prefix ending in value 0; f[2]=ending in value 1
  f <- c(0, 0); if (!is.na(x[1])) { if (x[1]==0) f[2] <- INF else f[1] <- INF }
  for (j in 2:L) {
    nf <- c(cmb(f[1] + 1, f[2] + 0),           # new value 0
            cmb(f[2] + 1, f[1] + 0))           # new value 1
    if (!is.na(x[j])) { if (x[j]==0) nf[2] <- INF else nf[1] <- INF }
    f <- nf
  }
  as.integer(cmb(f[1], f[2]))
}

# ---- (b) brute force over all 2^(#NA) fillings ----
chain_brute <- function(x, want = "max") {
  L <- length(x); if (L <= 1) return(0L)
  na <- which(is.na(x)); k <- length(na)
  agr <- function(v) sum(v[-L] == v[-1])
  if (k == 0) return(agr(x))
  best <- if (want=="max") -Inf else Inf; cmb <- if (want=="max") max else min
  for (m in 0:(2^k - 1)) {
    v <- x; bits <- as.integer(intToBits(m))[1:k]; v[na] <- bits
    best <- cmb(best, agr(v))
  }
  as.integer(best)
}

# ---- (c) VERIFY DP == brute force on many random small chains ----
cat("=== Verify single-chain DP vs brute force ===\n")
fails <- 0; ntest <- 4000
for (i in 1:ntest) {
  L <- sample(1:9, 1)
  x <- sample(c(0L,1L,NA), L, replace = TRUE, prob = c(.4,.4,.2))
  for (w in c("max","min")) {
    if (chain_opt(x, w) != chain_brute(x, w)) {
      fails <- fails + 1
      if (fails <= 5) cat("  MISMATCH", w, ": x =", paste(ifelse(is.na(x),".",x),collapse=""),
                          " dp=", chain_opt(x,w), " brute=", chain_brute(x,w), "\n")
    }
  }
}
cat(sprintf("single-chain: %d/%d mismatches\n\n", fails, ntest*2))

# ---- (d) vectorised multi-column DP (all M columns of a day-matrix at once) ----
# X: D x M matrix (rows = consecutive days), entries 0/1/NA.
# adj: logical length D-1, TRUE if days d,d+1 are consecutive (chain connected).
# returns total max (or min) agreements summed over all columns & connected pairs.
matrix_opt <- function(X, adj, want = "max") {
  D <- nrow(X); M <- ncol(X)
  INF <- if (want=="max") -Inf else Inf
  pick <- if (want=="max") pmax else pmin
  init_row <- function(r) {
    f0 <- rep(0, M); f1 <- rep(0, M)
    f0[which(r==1)] <- INF; f1[which(r==0)] <- INF
    list(f0=f0, f1=f1)
  }
  st <- init_row(X[1,]); total <- rep(0, M)
  for (d in 2:D) {
    if (!adj[d-1]) {                              # chain breaks: bank this run, restart
      total <- total + pick(st$f0, st$f1)
      st <- init_row(X[d,]); next
    }
    nf0 <- pick(st$f0 + 1, st$f1 + 0)
    nf1 <- pick(st$f1 + 1, st$f0 + 0)
    r <- X[d,]; nf0[which(r==1)] <- INF; nf1[which(r==0)] <- INF
    st <- list(f0=nf0, f1=nf1)
  }
  total <- total + pick(st$f0, st$f1)
  sum(total)
}

# ---- (e) VERIFY matrix_opt vs summing chain_brute column-by-column ----
cat("=== Verify vectorised matrix DP vs per-column brute force ===\n")
mfails <- 0; mtest <- 600
for (i in 1:mtest) {
  D <- sample(2:6,1); M <- sample(1:5,1)
  X <- matrix(sample(c(0L,1L,NA), D*M, replace=TRUE, prob=c(.4,.4,.2)), D, M)
  adj <- sample(c(TRUE,FALSE), D-1, replace=TRUE, prob=c(.7,.3))
  for (w in c("max","min")) {
    # brute: split into connected runs, brute each column within each run
    runs <- split(1:D, cumsum(c(0, !adj)))
    bt <- 0
    for (t in 1:M) for (rr in runs) bt <- bt + chain_brute(X[rr, t], w)
    dp <- matrix_opt(X, adj, w)
    if (dp != bt) { mfails <- mfails+1
      if (mfails<=5) cat("  MISMATCH", w, " dp=",dp," brute=",bt,"\n") }
  }
}
cat(sprintf("matrix DP: %d/%d mismatches\n\n", mfails, mtest*2))

if (fails==0 && mfails==0) cat("ALL VERIFICATIONS PASSED — exact bounds algorithm is correct.\n") else cat("!!! FAILURES PRESENT — fix needed.\n")
