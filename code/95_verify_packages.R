# =====================================================================
# VERIFY that GGIR and sleepreg SRI are computed correctly, by feeding
# synthetic participants whose SRI is known analytically, and by cross-
# checking against the hand formula on the SAME series.
#   Case A  7 identical days                      -> SRI = +100
#   Case B  7 days alternating pattern/complement -> SRI = -100
#   Case C  consecutive days agree on half epochs -> SRI =   0
#   Case D  random-ish fixed pattern              -> hand == packages
# =====================================================================
suppressMessages({library(data.table); library(GGIR)})
options(width=150)
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())
rdir  <- file.path(base,"01_recompute_SRI")
source(file.path(rdir,"sleepreg_src","sleepreg_allFunctions.R")); load(file.path(rdir,"sleepreg_src","sysdata.rda"))

M <- 1440L; D <- 7L
# build a participant from a D x 1440 matrix of 0/1 (rows=days). timestamp = real clock, midnight-aligned.
build <- function(mat){
  s  <- as.integer(t(mat))                       # row-major: day1 min1..1440, day2...
  ts <- 43200 + (0:(length(s)-1))*60                      # seconds from 1970-01-01 00:00:00, 60s steps
  data.table(binary=s, ts=ts,
             day=rep(1:nrow(mat), each=ncol(mat)),
             dow=rep(((0:(nrow(mat)-1))%%7)+1, each=ncol(mat)))   # day-of-week 1..7
}
hand_sri <- function(s,d,M=1440L){ n<-length(s); i<-1:(n-M); cur<-d[i]; nx<-d[i+M]
  con<-(nx-cur==1)|(cur==7&nx==1); a<-s[i]; b<-s[i+M]; v<-con&!is.na(a)&!is.na(b)
  Nv<-sum(v); -100+200/Nv*sum(a[v]==b[v]) }
ggir_sri <- function(ts,s){ df<-data.frame(time=as.POSIXct(ts,origin="1970-01-01",tz="UTC"),
  sleepstate=as.numeric(s), invalid=as.integer(is.na(s)))
  o<-GGIR::CalcSleepRegularityIndex(data=df,epochsize=60,desiredtz="UTC"); mean(o$SleepRegularityIndex,na.rm=TRUE) }
sleepreg_sri <- function(dt){ td<-file.path(tempdir(),paste0("v",sample(1e6,1))); dir.create(td)
  fwrite(dt[,.(binary,ts)], file.path(td,"p.csv"))
  SWS_from_binarySW(binarySWdir=td, tsCol=2, binaryCol=1)
  invisible(capture.output(SRI_from_binary(binarydir=file.path(td,"SWS_output"),
    alloutdir=file.path(td,"out"), col.trans=1,col.timestamp=2, wr.raster=FALSE, overwr=TRUE, minSRIdays=3)))
  as.numeric(fread(file.path(td,"out","SRI.csv"))$SRI[1]) }

run <- function(name, mat, expected){
  dt <- build(mat)
  h <- hand_sri(dt$binary, dt$dow)
  g <- ggir_sri(dt$ts, dt$binary)
  s <- sleepreg_sri(dt)
  cat(sprintf("%-42s expected=%6.1f | hand=%7.2f  GGIR=%7.2f  sleepreg=%7.2f\n", name, expected, h, g, s))
}

# patterns (each day = 1440 epochs; 1=sleep,0=wake)
day_pat <- c(rep(1L,480), rep(0L,960))             # sleep 8h then wake 16h
# A: identical every day
A <- matrix(rep(day_pat, D), nrow=D, byrow=TRUE)
# B: alternate pattern / complement
B <- do.call(rbind, lapply(0:(D-1), function(k) if(k%%2==0) day_pat else 1L-day_pat))
# C: agree on second half (720-1439), disagree on first half across consecutive days
half_even <- c(rep(1L,720), rep(0L,720)); half_odd <- c(rep(0L,720), rep(0L,720))
C <- do.call(rbind, lapply(0:(D-1), function(k) if(k%%2==0) half_even else half_odd))
# D: fixed pseudo-random pattern, same each-day-shifted; compute hand as reference
set.seed(7); base_pat <- rbinom(M,1,0.4)
Dm <- do.call(rbind, lapply(0:(D-1), function(k){ shift<-k*30; c(base_pat[(shift+1):M], base_pat[seq_len(shift)]) }))

cat("=== synthetic SRI verification (7 consecutive days, no missing) ===\n")
run("A: 7 identical days", A, 100)
run("B: alternating pattern/complement", B, -100)
run("C: agree on half the epochs", C, 0)
dt <- build(Dm); cat(sprintf("%-42s hand=%7.2f  GGIR=%7.2f  sleepreg=%7.2f  (hand is reference)\n",
   "D: shifted pseudo-random", hand_sri(dt$binary,dt$dow), ggir_sri(dt$ts,dt$binary), sleepreg_sri(dt)))
cat("\nPASS if GGIR & sleepreg match expected (A/B/C) and match hand (D).\n")
