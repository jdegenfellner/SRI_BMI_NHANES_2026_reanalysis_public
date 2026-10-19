# make_crosscheck_list.py [n_je_zyklus] -> aws_list_crosscheck.tsv : bereits lokal fertige Teilnehmer fuer den Plattformvergleich (Mac vs. Linux)
import json, os, sys
n = int(sys.argv[1]) if len(sys.argv) > 1 else 10; base = "/Volumes/LaCie/NHANES_1112_1314"
L = json.load(open(f"{base}/_download/remote_listing.json")); got = {"pax_g": 0, "pax_h": 0}; out = []
for l in open(f"{base}/_work/analytic_random_order.tsv").read().splitlines()[1:]:
    r, s, c = l.split("\t")
    if got[c] < n and os.path.exists(f"{base}/_ggir/{c}/done/{s}.txt") and os.path.exists(f"{base}/_swan/{c}/done/{s}.txt"): out.append((r, s, c, L[s]["size"])); got[c] += 1
open("aws_list_crosscheck.tsv", "w").write("rank\tseqn\tcycle\tsize\n" + "".join("\t".join(map(str, o)) + "\n" for o in out)); print(len(out), "Teilnehmer")
