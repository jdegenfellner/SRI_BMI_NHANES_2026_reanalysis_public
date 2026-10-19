#!/bin/bash
# one_participant.sh <seqn> <pax_g|pax_h> <size>  (im Container)
# Archiv von der CDC laden, GGIR und SWaN parallel rechnen (wie lokal), Ergebnisse nach $NH_BASE/_ggir und $NH_BASE/_swan, Archiv loeschen.
seqn=$1; cyc=$2; size=$3; export NH_BASE=${NH_BASE:-/work/out}; export NH_TMP=${NH_TMP:-/work/tmp}
raw=/work/raw/$seqn.tar.bz2; mkdir -p /work/raw "$NH_BASE/_ggir/$cyc" "$NH_BASE/_swan/$cyc" "$NH_TMP"; log=$NH_BASE/cloud.log
[ -f "$NH_BASE/_ggir/$cyc/done/$seqn.txt" ] && [ -f "$NH_BASE/_swan/$cyc/done/$seqn.txt" ] && exit 0
[ -f /work/STOP ] && exit 0
t0=$(date +%s)
for try in 1 2 3 4 5 6; do
  curl -s -A "Mozilla/5.0" --retry 3 --retry-delay 10 --connect-timeout 30 --speed-time 120 --speed-limit 1000 -C - -o "$raw.part" "https://ftp.cdc.gov/pub/$cyc/$seqn.tar.bz2"
  [ "$(stat -c %s "$raw.part" 2>/dev/null)" = "$size" ] && { mv "$raw.part" "$raw"; break; }
  [ "$(stat -c %s "$raw.part" 2>/dev/null || echo 0)" -gt "$size" ] && rm -f "$raw.part"; sleep $((try*20))
done
[ -f "$raw" ] || { echo "$(date +%FT%T) $seqn $cyc DOWNLOAD FAIL" >> "$log"; exit 1; }
t1=$(date +%s); export NH_RAW_FILE=$raw
if [ ! -f "$NH_BASE/_ggir/$cyc/done/$seqn.txt" ]; then Rscript /opt/pipeline/ggir_one_cloud.R "$seqn" "$cyc" > "/work/tmp/ggir_$seqn.log" 2>&1 & gp=$!; fi
if [ ! -f "$NH_BASE/_swan/$cyc/done/$seqn.txt" ]; then /opt/pipeline/swan_one_cloud.sh "$seqn" "$cyc" > "/work/tmp/swan_$seqn.log" 2>&1 & sp=$!; fi
[ -n "$gp" ] && wait $gp; [ -n "$sp" ] && wait $sp
t2=$(date +%s); ok=ok; [ -f "$NH_BASE/_ggir/$cyc/done/$seqn.txt" ] && [ -f "$NH_BASE/_swan/$cyc/done/$seqn.txt" ] || ok=FEHLT
echo "$(date +%FT%T) $seqn $cyc $ok download $((t1-t0)) s, rechnen $((t2-t1)) s" >> "$log"
[ "$ok" = ok ] && rm -f "$raw" "/work/tmp/ggir_$seqn.log" "/work/tmp/swan_$seqn.log"
