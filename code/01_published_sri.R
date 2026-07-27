# =====================================================================
# 01_published_sri.R
# The baseline for this re-analysis is the SRI exactly as published in the
# original study (Degenfellner et al. 2026, J Biol Rhythms). Those values were
# computed once and frozen as dated snapshots (Aug 2024) and used throughout the
# publication. The two snapshots are shipped here so this re-analysis uses the
# EXACT published values.
#
# Note on reproducibility: re-running the original preparation code on the
# current NHANES data does NOT reproduce these values bit-for-bit (the original
# code base drifted after the snapshot was frozen). The shipped snapshot is
# therefore the authoritative source of the published SRI.
#
# Output: data/published_sri.csv  (SEQN, SRI_published)
# =====================================================================
suppressMessages(library(data.table))
base <- Sys.getenv("SRI_PROJECT_ROOT", getwd())
source(file.path(base, if (dir.exists(file.path(base, "code"))) "code" else "01_recompute_SRI", "00_paths.R"))

s12 <- as.data.table(readRDS(insrc("published_SRI_2011_12.RDS")))
s34 <- as.data.table(readRDS(insrc("published_SRI_2013_14.RDS")))
published <- rbind(s12[, .(SEQN, SRI_published = SRI_no_imputation)],
                   s34[, .(SEQN, SRI_published = SRI_no_imputation)])
# canonical, minimal baseline file (only the relevant columns), both formats
fwrite(published,  file.path(inter, "published_sri.csv"))
saveRDS(published, file.path(inter, "published_sri.rds"))
cat("published SRI:", nrow(published), "participants -> published_sri.{csv,rds}\n")

# Verification: the shipped snapshot equals the "Original / published" baseline
# used everywhere in this re-analysis.
chk <- as.data.table(readRDS(insrc("sri_all_methods.rds")))[!is.na(SRI_orig)]
m   <- merge(chk[, .(SEQN, SRI_orig)], published, by = "SEQN")
cat(sprintf("verification: snapshot == analysis baseline for %d/%d participants (cor = %.4f)\n",
            sum(abs(m$SRI_orig - m$SRI_published) < 1e-6), nrow(m),
            cor(m$SRI_orig, m$SRI_published)))
