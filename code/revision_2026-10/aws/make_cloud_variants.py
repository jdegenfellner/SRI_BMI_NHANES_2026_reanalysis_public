# Erzeugt ggir_one_cloud.R und swan_one_cloud.sh aus den lokal verwendeten Skripten.
# Nur Pfade und Hilfsprogramme werden ersetzt; die GGIR()-Parameter und der SWaN-Python-Block bleiben Zeichen fuer Zeichen gleich.
import sys, difflib
W = sys.argv[1]; out = sys.argv[2]
def patch(src, dst, reps):
    s = open(f"{W}/{src}").read()
    for a, b in reps:
        assert s.count(a) == 1, (src, a, s.count(a))
        s = s.replace(a, b)
    open(f"{out}/{dst}", "w").write(s)
    print(f"== {dst}: " + "".join(l for l in difflib.unified_diff(open(f"{W}/{src}").read().splitlines(1), s.splitlines(1), n=0) if l.startswith(("+", "-")) and not l.startswith(("+++", "---"))))
patch("ggir_one.R", "ggir_one_cloud.R", [
  ('base <- "/Volumes/LaCie/NHANES_1112_1314"', 'base <- Sys.getenv("NH_BASE", "/Volumes/LaCie/NHANES_1112_1314")'),
  ('raw <- file.path(base, if (cyc == "pax_g") "raw_accel_data_201112" else "raw_accel_data_201314", paste0(seqn, ".tar.bz2"))',
   'raw <- Sys.getenv("NH_RAW_FILE", file.path(base, if (cyc == "pax_g") "raw_accel_data_201112" else "raw_accel_data_201314", paste0(seqn, ".tar.bz2")))'),
  ('W <- file.path(Sys.getenv("HOME"), ".ggir_tmp", paste0("ggir_", seqn))',
   'W <- file.path(Sys.getenv("NH_TMP", file.path(Sys.getenv("HOME"), ".ggir_tmp")), paste0("ggir_", seqn))'),
  ('t1 <- Sys.time()\nGGIR(mode = 1:5,', 't1 <- Sys.time()\nif (nzchar(Sys.getenv("NH_GGIR_PATCH"))) source(Sys.getenv("NH_GGIR_PATCH"))\nGGIR(mode = 1:5,')])
patch("swan_one.sh", "swan_one_cloud.sh", [
  ('base=/Volumes/LaCie/NHANES_1112_1314; PY=$HOME/.venvs/swan38/bin/python', 'base=${NH_BASE:-/Volumes/LaCie/NHANES_1112_1314}; PY=${NH_PY:-$HOME/.venvs/swan38/bin/python}'),
  ('[ "$cyc" = "pax_h" ] && raw=$base/raw_accel_data_201314/$seqn.tar.bz2', '[ "$cyc" = "pax_h" ] && raw=$base/raw_accel_data_201314/$seqn.tar.bz2; raw=${NH_RAW_FILE:-$raw}'),
  ('T=$HOME/.ggir_tmp/swan_$seqn', 'T=${NH_TMP:-$HOME/.ggir_tmp}/swan_$seqn'),
  ('"$base/_work/swan_prep.R"', '"${NH_SWAN_PREP:-$base/_work/swan_prep.R}"'),
  ('gzcat "$T/final.csv.gz"', '${NH_ZCAT:-gzcat} "$T/final.csv.gz"')])
