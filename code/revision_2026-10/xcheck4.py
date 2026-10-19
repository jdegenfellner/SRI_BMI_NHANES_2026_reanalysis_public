import importlib.util, sys, numpy as np, pandas as pd
# usage: python3 xcheck4.py <path to pyActigraphy/sleep/scoring/sri.py>  (run in the folder with xcheck4_series.csv)
spec = importlib.util.spec_from_file_location("sri", sys.argv[1] if len(sys.argv) > 1 else "pyact312/lib/python3.12/site-packages/pyActigraphy/sleep/scoring/sri.py"); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
df = pd.read_csv("xcheck4_series.csv"); rows = []
for seqn, g in df.groupby("SEQN", sort=True):
    g = g.sort_values("ts"); ts = pd.Series(g["s"].values.astype(float), index=pd.to_datetime(g["ts"].values, unit="s")); ts.index.freq = "1min"; rows.append((seqn, m.sri(ts)))
pd.DataFrame(rows, columns=["SEQN", "SRI_pyact"]).to_csv("xcheck4_pyact.csv", index=False); print("py done", len(rows))
