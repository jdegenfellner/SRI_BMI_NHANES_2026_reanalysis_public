#!/usr/bin/env python3
# ingest_staging.py <staging_run_dir> : uebernimmt fertige Cloud-Teilnehmer aus der lokalen Sicherung (LaCie/_cloud/<Instanz>/<Lauf>) in den Hauptbestand.
# Nur Teilnehmer mit GGIR- UND SWaN-Marker, die im Hauptbestand noch nicht fertig sind. Ergebnisdateien zuerst, Done-Marker zuletzt
# (Dispatcher und Auswertung sehen den Teilnehmer erst danach als fertig). Protokoll: _work/cloud_ingested.log (Plattformherkunft).
import subprocess, sys, os, datetime
src = sys.argv[1].rstrip("/"); base = "/Volumes/LaCie/NHANES_1112_1314"; tag = "/".join(src.split("/")[-2:])
n = 0
for cyc in ("pax_g", "pax_h"):
    ls = lambda p: {f[:-4] for f in os.listdir(p) if f.endswith(".txt") and not f.startswith(".")} if os.path.isdir(p) else set()
    for seqn in sorted(ls(f"{src}/_ggir/{cyc}/done") & ls(f"{src}/_swan/{cyc}/done")):
        if os.path.exists(f"{base}/_ggir/{cyc}/done/{seqn}.txt") and os.path.exists(f"{base}/_swan/{cyc}/done/{seqn}.txt"): continue
        for t in ("_ggir", "_swan"): os.makedirs(f"{base}/{t}/{cyc}/done", exist_ok=True)
        subprocess.run(["rsync", "-a", f"{src}/_ggir/{cyc}/{seqn}/", f"{base}/_ggir/{cyc}/{seqn}/"], check=True)
        subprocess.run(["rsync", "-a", f"{src}/_swan/{cyc}/{seqn}_final.csv.gz", f"{base}/_swan/{cyc}/"], check=True)
        subprocess.run(["rsync", "-a", f"{src}/_swan/{cyc}/done/{seqn}.txt", f"{base}/_swan/{cyc}/done/"], check=True)
        subprocess.run(["rsync", "-a", f"{src}/_ggir/{cyc}/done/{seqn}.txt", f"{base}/_ggir/{cyc}/done/"], check=True)
        open(f"{base}/_work/cloud_ingested.log", "a").write(f"{datetime.datetime.now():%FT%T} {seqn} {cyc} {tag}\n"); n += 1
print(n, "Teilnehmer uebernommen")
