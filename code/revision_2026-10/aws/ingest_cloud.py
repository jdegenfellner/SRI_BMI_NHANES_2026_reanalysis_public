#!/usr/bin/env python3
# ingest_cloud.py <ip> [remote_base=/work/prod] : uebernimmt fertige Cloud-Teilnehmer auf die LaCie.
# Nur Teilnehmer mit GGIR- UND SWaN-Marker in der Cloud, die lokal noch nicht fertig sind. Ergebnisdateien zuerst, Done-Marker zuletzt
# (der lokale Dispatcher und die Auswertung sehen den Teilnehmer erst danach als fertig). Protokoll: _work/cloud_ingested.log (Plattformherkunft).
import subprocess, sys, os, datetime
ip = sys.argv[1]; rb = sys.argv[2] if len(sys.argv) > 2 else "/work/prod"; base = os.environ.get("RAW_ROOT", "/Volumes/LaCie/NHANES_1112_1314")
ssh = ["ssh", "-i", os.path.expanduser(os.environ.get("CLOUD_SSH_KEY", "~/.ssh/cloud-key.pem")), "-o", "StrictHostKeyChecking=accept-new", "-o", "BatchMode=yes"]
rs = ["rsync", "-a", "-e", " ".join(ssh)]
def remote(cmd): return subprocess.run(ssh + [f"ubuntu@{ip}", cmd], capture_output=True, text=True, check=True).stdout.split()
def run(a): subprocess.run(a, check=True)
n = 0
for cyc in ("pax_g", "pax_h"):
    g = {f[:-4] for f in remote(f"ls {rb}/_ggir/{cyc}/done 2>/dev/null")}; s = {f[:-4] for f in remote(f"ls {rb}/_swan/{cyc}/done 2>/dev/null")}
    for seqn in sorted(g & s):
        if os.path.exists(f"{base}/_ggir/{cyc}/done/{seqn}.txt") and os.path.exists(f"{base}/_swan/{cyc}/done/{seqn}.txt"): continue
        for t in ("_ggir", "_swan"): os.makedirs(f"{base}/{t}/{cyc}/done", exist_ok=True)
        run(rs + [f"ubuntu@{ip}:{rb}/_ggir/{cyc}/{seqn}/", f"{base}/_ggir/{cyc}/{seqn}/"])
        run(rs + [f"ubuntu@{ip}:{rb}/_swan/{cyc}/{seqn}_final.csv.gz", f"{base}/_swan/{cyc}/"])
        run(rs + [f"ubuntu@{ip}:{rb}/_swan/{cyc}/done/{seqn}.txt", f"{base}/_swan/{cyc}/done/"])
        run(rs + [f"ubuntu@{ip}:{rb}/_ggir/{cyc}/done/{seqn}.txt", f"{base}/_ggir/{cyc}/done/"])
        open(f"{base}/_work/cloud_ingested.log", "a").write(f"{datetime.datetime.now():%FT%T} {seqn} {cyc} {ip}\n"); n += 1
print(n, "Teilnehmer uebernommen")
