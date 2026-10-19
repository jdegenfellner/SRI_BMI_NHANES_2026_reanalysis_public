suppressMessages({library(data.table); library(haven)})
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())   # working tree root (contains 00_original_paper_materials/ and 01_recompute_SRI/)
hdr <- as.data.table(read_xpt(file.path(base,"01_recompute_SRI/data_local/PAXHD_G_2011_12.XPT")))
cat("PAXHD_G columns:", paste(names(hdr), collapse=", "), "\n")
print(head(hdr, 4))
cat("class PAXFTIME:", class(hdr$PAXFTIME), " unique first 10:", paste(head(unique(hdr$PAXFTIME),10), collapse=" | "), "\n\n")
an <- unique(as.data.table(readRDS(file.path(base,"00_original_paper_materials/data/df_an1_imputed_MAIN.RDS")))$SEQN)
ser <- as.data.table(readRDS(file.path(base,"01_recompute_SRI/intermediate/series_combined.rds")))[SEQN %in% an]
setorder(ser, SEQN, UNIX_TIMESTAMP)
# TRUE clock time: every retained day has exactly 1440 minutes and is a calendar day, so its first minute is 00:00.
ser[, day_start := UNIX_TIMESTAMP[1], by = .(SEQN, PAXDAYM)]
ser[, tod_h := (UNIX_TIMESTAMP - day_start)/3600]
stopifnot(all(ser[, .N, by=.(SEQN,PAXDAYM)]$N == 1440))
cat("range tod_h:", range(ser$tod_h), "\n")
ser[, state := fifelse(is.na(binary_sleep), -1, binary_sleep)]
ser[, hour := floor(tod_h)]
byh <- ser[state != -1, .(p_sleep = round(mean(state == 1),3)), by = hour][order(hour)]
cat("\n== share classified asleep by TRUE clock hour ==\n"); print(byh)
cat("share of all sleep minutes 10:00-20:00:", round(ser[state==1, mean(tod_h >= 10 & tod_h < 20)],3), "\n")
cat("share of all sleep minutes 22:00-08:00:", round(ser[state==1, mean(tod_h >= 22 | tod_h < 8)],3), "\n")
# night-window fragmentation 23:00-06:00 (within calendar day, so 23-24 and 0-6 belong to different PAXDAYM; use noon-to-noon nights)
ser[, night_id := fifelse(tod_h >= 12, as.integer(PAXDAYM), as.integer(PAXDAYM) - 1L)]
night <- ser[state != -1 & (tod_h >= 23 | tod_h < 6)]
nt <- night[, { s <- state; .(n_trans = sum(diff(s) != 0), p_sleep = mean(s==1), n = .N) }, by = .(SEQN, night_id)][n >= 300]
cat("\n== 23:00-06:00 window: sleep<->wake transitions per night ==\n"); print(round(quantile(nt$n_trans, c(.1,.25,.5,.75,.9,.95)),1))
cat("mean p(sleep) 23-06:", round(mean(nt$p_sleep),3), " median:", round(median(nt$p_sleep),3), "\n")
# wake bouts inside 23-06 window: lengths
nr <- night[, { r <- rle(state); .(st = r$values, len = r$lengths) }, by = .(SEQN, night_id)]
wn <- nr[st == 0]
cat("== wake bout lengths inside 23-06 (min) ==\n"); print(quantile(wn$len, c(.25,.5,.75,.9,.95)))
cat("share of night wake MINUTES in bouts <=5 min:", round(sum(wn$len[wn$len<=5])/sum(wn$len),3), " <=15:", round(sum(wn$len[wn$len<=15])/sum(wn$len),3), " >=30:", round(sum(wn$len[wn$len>=30])/sum(wn$len),3), "\n")
# daytime sleep: 10:00-20:00 sleep bouts
day <- ser[state != -1 & tod_h >= 10 & tod_h < 20]
dr <- day[, { r <- rle(state); .(st = r$values, len = r$lengths) }, by = .(SEQN, PAXDAYM)]
ds <- dr[st == 1]
cat("\n== daytime (10-20h) sleep bouts: per participant-day count and length ==\n")
cat("bouts per day median:", median(ds[, .N, by=.(SEQN,PAXDAYM)]$N), "\n"); print(quantile(ds$len, c(.25,.5,.75,.9,.95)))
cat("daytime sleep minutes per day, median:", median(ds[, sum(len), by=.(SEQN,PAXDAYM)]$V1), "\n")
saveRDS(byh, "byh_true.rds")
cat("DONE\n")
