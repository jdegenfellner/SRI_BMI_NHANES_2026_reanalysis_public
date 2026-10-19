suppressMessages(library(data.table))
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())   # working tree root (contains 00_original_paper_materials/ and 01_recompute_SRI/)
an <- unique(as.data.table(readRDS(file.path(base,"00_original_paper_materials/data/df_an1_imputed_MAIN.RDS")))$SEQN)
ser <- as.data.table(readRDS(file.path(base,"01_recompute_SRI/intermediate/series_combined.rds")))[SEQN %in% an]
setorder(ser, SEQN, UNIX_TIMESTAMP)
cat("participants:", uniqueN(ser$SEQN), " rows:", nrow(ser), "\n")
cat("PAXPREDM categories (minute shares):\n"); print(round(prop.table(table(ser$PAXPREDM)),4))
# time of day
ser[, tod_h := (UNIX_TIMESTAMP %% 86400)/3600]
# --- per participant-day: runs of binary_sleep (NA treated as its own state, then dropped) ---
ser[, state := fifelse(is.na(binary_sleep), -1, binary_sleep)]
runs <- ser[, {
  r <- rle(state)
  .(state = r$values, len = r$lengths)
}, by = .(SEQN, PAXDAYM)]
# transitions between sleep and wake per day (ignoring NA runs): count adjacent sleep<->wake pairs
trans <- ser[, {
  s <- state; s <- s[s != -1]
  .(n_sw_trans = sum(diff(s) != 0), n_sleep_min = sum(s == 1), n_wake_min = sum(s == 0), n_na = sum(state == -1))
}, by = .(SEQN, PAXDAYM)]
cat("\n== sleep<->wake transitions per participant-day (excluding NA minutes) ==\n")
print(round(quantile(trans$n_sw_trans, c(.05,.1,.25,.5,.75,.9,.95,.99)),1))
cat("mean:", round(mean(trans$n_sw_trans),1), "\n")
cat("sleep minutes/day median:", median(trans$n_sleep_min), " mean:", round(mean(trans$n_sleep_min),1), "\n")
# per participant mean transitions/day
pt <- trans[, .(trans_day = mean(n_sw_trans), sleep_min_day = mean(n_sleep_min)), by = SEQN]
cat("\n== per-participant mean transitions/day ==\n"); print(round(quantile(pt$trans_day, c(.05,.25,.5,.75,.95)),1))
# --- bout lengths ---
sb <- runs[state == 1]; wb <- runs[state == 0]
cat("\n== SLEEP bout length (min) ==\n"); print(quantile(sb$len, c(.1,.25,.5,.75,.9,.95,.99)))
cat("share of sleep MINUTES in bouts <=5 min:", round(sum(sb$len[sb$len<=5])/sum(sb$len),3),
    " <=15:", round(sum(sb$len[sb$len<=15])/sum(sb$len),3),
    " <=30:", round(sum(sb$len[sb$len<=30])/sum(sb$len),3),
    " >=120:", round(sum(sb$len[sb$len>=120])/sum(sb$len),3), "\n")
cat("share of sleep BOUTS <=5 min:", round(mean(sb$len<=5),3), " <=15:", round(mean(sb$len<=15),3), "\n")
cat("\n== WAKE bout length (min) ==\n"); print(quantile(wb$len, c(.1,.25,.5,.75,.9,.95,.99)))
cat("share of wake MINUTES in bouts <=5 min:", round(sum(wb$len[wb$len<=5])/sum(wb$len),3),
    " <=15:", round(sum(wb$len[wb$len<=15])/sum(wb$len),3), "\n")
# --- sleep by time of day (daytime sleep = misclassification or naps?) ---
ser[, hour := floor(tod_h)]
byh <- ser[state != -1, .(p_sleep = mean(state == 1)), by = hour][order(hour)]
cat("\n== share classified asleep by clock hour ==\n"); print(byh[, .(hour, p_sleep = round(p_sleep,3))])
cat("share of all sleep minutes falling 10:00-20:00:", round(ser[state==1, mean(tod_h >= 10 & tod_h < 20)],3), "\n")
# --- longest sleep bout per day (main sleep period) and sleep outside it ---
lsb <- runs[state == 1, .(longest = max(len), n_bouts = .N, total_sleep = sum(len)), by = .(SEQN, PAXDAYM)]
cat("\n== longest sleep bout per day (min) ==\n"); print(quantile(lsb$longest, c(.1,.25,.5,.75,.9)))
cat("n sleep bouts/day median:", median(lsb$n_bouts), " mean:", round(mean(lsb$n_bouts),1), "\n")
cat("share of sleep minutes outside the longest bout:", round(1 - sum(lsb$longest)/sum(lsb$total_sleep),3), "\n")
# --- wake bouts INSIDE the night (between 23:00 and 06:00): WASO-like fragmentation ---
night <- ser[state != -1 & (tod_h >= 23 | tod_h < 6)]
nt <- night[, { s <- state; .(n_trans = sum(diff(s) != 0), p_sleep = mean(s==1)) }, by = .(SEQN, PAXDAYM)]
cat("\n== 23:00-06:00 window: transitions per night ==\n"); print(round(quantile(nt$n_trans, c(.1,.25,.5,.75,.9,.95)),1))
cat("mean p(sleep) 23-06:", round(mean(nt$p_sleep),3), "\n")
saveRDS(list(trans=trans, pt=pt, sb=sb, wb=wb, byh=byh, lsb=lsb, nt=nt), "frag_summary.rds")
cat("\nDONE\n")
