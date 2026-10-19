# make_list.py <rank_von> <rank_bis> -> aws_list_<von>_<bis>.tsv (rank, seqn, cycle, size) fuer die Teilnehmer, die lokal noch nicht fertig sind
import json, sys, os
von, bis = int(sys.argv[1]), int(sys.argv[2]); base = "/Volumes/LaCie/NHANES_1112_1314"
L = json.load(open(f"{base}/_download/remote_listing.json"))
done = lambda c, s: os.path.exists(f"{base}/_ggir/{c}/done/{s}.txt") and os.path.exists(f"{base}/_swan/{c}/done/{s}.txt")
excl = set(open(f"{base}/_work/excluded_seqn.txt").read().split())
rows = [l.split("\t") for l in open(f"{base}/_work/analytic_random_order.tsv").read().splitlines()[1:]]
out = [(int(r[0]), r[1], r[2], L[r[1]]["size"]) for r in rows if von <= int(r[0]) <= bis and r[1] not in excl and not done(r[2], r[1])]
fn = f"aws_list_{von}_{bis}.tsv"; open(fn, "w").write("rank\tseqn\tcycle\tsize\n" + "".join("\t".join(map(str, o)) + "\n" for o in out)); print(fn, len(out), "Teilnehmer,", round(sum(o[3] for o in out) / 1e9, 1), "GB")
