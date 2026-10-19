# Sample selection shared by raw_derivation_summary.R and raw_derivation_models.R
# pick_sample(res, ord, mode, rmax, excluded)
#   res       participant-level results (raw_derivation_participants.rds)
#   ord       data.table(SEQN, rank) : fixed random order of all analytic participants (analytic_random_order.tsv)
#   mode      "prefix": participants with rank <= R and all four series valid, R = largest rank such that every participant with a smaller
#             or equal rank has all four series processed or is listed in `excluded` (a random sample of the analytic sample);
#             "all": every participant with all four series valid (final analysis)
pick_sample <- function(res, ord, mode = "prefix", rmax = Inf, excluded = integer()) {
  cnt <- res[, .(n_deriv = uniqueN(deriv), ok = all(!is.na(sri_A))), by = SEQN]
  complete <- cnt[n_deriv == 4 & ok, SEQN]
  ord <- copy(ord); ord[, complete := SEQN %in% complete]; ord[, resolved := SEQN %in% cnt[n_deriv == 4, SEQN] | SEQN %in% excluded]
  r_avail <- if (all(ord$resolved)) nrow(ord) else min(ord[resolved == FALSE, rank]) - 1L
  if (mode == "prefix") { R <- min(r_avail, rmax, nrow(ord)); samp <- ord[rank <= R & complete, SEQN] } else { R <- NA_integer_; samp <- complete }
  list(samp = samp, R = R, complete = complete, ord = ord, n_failed = if (mode == "prefix") sum(ord[rank <= R, SEQN] %in% excluded) else sum(ord$SEQN %in% excluded))
}
