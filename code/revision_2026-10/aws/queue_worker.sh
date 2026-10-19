#!/bin/bash
# queue_worker.sh <liste.tsv> <slots>  (im Container). Liste: rank<TAB>seqn<TAB>cycle<TAB>size, in Zufallsreihenfolge. Fortsetzbar (fertige werden uebersprungen).
list=$1; slots=$2
grep -v '^rank' "$list" | awk -F'\t' '{print $2, $3, $4}' | xargs -P "$slots" -L 1 /opt/pipeline/one_participant.sh
echo "$(date +%FT%T) QUEUE FERTIG" >> ${NH_BASE:-/work/out}/cloud.log
