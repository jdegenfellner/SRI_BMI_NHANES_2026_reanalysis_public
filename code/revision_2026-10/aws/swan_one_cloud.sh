#!/bin/bash
# swan_one.sh <seqn> <pax_g|pax_h> : SWaN (oeffentliches Paket 1.18, g, ohne ISMR) -> _swan/<cycle>/<seqn>_final.csv.gz
seqn=$1; cyc=$2; base=${NH_BASE:-/Volumes/LaCie/NHANES_1112_1314}; PY=${NH_PY:-$HOME/.venvs/swan38/bin/python}
out=$base/_swan/$cyc; mkdir -p "$out/done"; [ -f "$out/done/$seqn.txt" ] && { echo "$seqn schon fertig"; exit 0; }
raw=$base/raw_accel_data_201112/$seqn.tar.bz2; [ "$cyc" = "pax_h" ] && raw=$base/raw_accel_data_201314/$seqn.tar.bz2; raw=${NH_RAW_FILE:-$raw}
T=${NH_TMP:-$HOME/.ggir_tmp}/swan_$seqn; rm -rf "$T"; mkdir -p "$T/extract" "$T/mhealth"; t0=$(date +%s)
tar -xjf "$raw" -C "$T/extract" && nice -n 10 Rscript "${NH_SWAN_PREP:-$base/_work/swan_prep.R}" "$T/extract" "$T/mhealth" && rm -rf "$T/extract"; t1=$(date +%s)
nice -n 10 perl -e 'alarm shift; exec @ARGV' 3600 $PY - "$T/mhealth" <<'PYEOF'
import sys, os, glob, warnings, pandas as pd
warnings.simplefilter("ignore")
from SWaN_accel import swan_first_pass, swan_second_pass
root = sys.argv[1]
for f in sorted(glob.glob(os.path.join(root, "*", "*", "*.sensor.csv"))):
    df = pd.read_csv(f, parse_dates=["HEADER_TIME_STAMP"]); swan_first_pass.main(df=df, file_path=f.replace(".sensor.csv", ".featureOri0to180.csv.gz"), sampling_rate=80); os.remove(f)
days = sorted(d for d in os.listdir(root) if os.path.isdir(os.path.join(root, d)))
for d in days: swan_second_pass.main(day_folder=os.path.join(root, d))
fin = [pd.read_csv(os.path.join(root, d, f"SWaN_{d}_final.csv")) for d in days if os.path.exists(os.path.join(root, d, f"SWaN_{d}_final.csv"))]
pd.concat(fin, ignore_index=True).to_csv(os.path.join(root, "..", "final.csv.gz"), index=False, compression="gzip"); print("Tage", len(fin))
PYEOF
t2=$(date +%s); [ -f "$T/final.csv.gz" ] || { echo "$(date +%FT%T) $seqn $cyc SWaN ohne Ergebnis nach $((t2-t1)) s (Zeitlimit 3600 s oder Fehler)" >> $base/_swan/failed.log; rm -rf "$T"; exit 1; }; raw_start=$(${NH_ZCAT:-gzcat} "$T/final.csv.gz" | sed -n 2p | cut -d, -f1)
cp "$T/final.csv.gz" "$out/${seqn}_final.csv.gz" && printf "seqn\t%s\ncycle\t%s\nraw_start\t%s\nprep_s\t%s\nswan_s\t%s\n" "$seqn" "$cyc" "$raw_start" $((t1-t0)) $((t2-t1)) > "$out/done/$seqn.txt"
rm -rf "$T"; echo "$seqn fertig: prep $((t1-t0)) s, SWaN $((t2-t1)) s"
