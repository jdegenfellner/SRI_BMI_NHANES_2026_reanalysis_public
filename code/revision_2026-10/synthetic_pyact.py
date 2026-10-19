# pyActigraphy sri() on the synthetic recordings (complete seven days), compared with the Phillips formula.
import sys, importlib.util, numpy as np, pandas as pd
venv_sri, cases_csv, out_csv = sys.argv[1], sys.argv[2], sys.argv[3]
spec = importlib.util.spec_from_file_location("sri", venv_sri); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
df = pd.read_csv(cases_csv); rows = []
for case, g in df.groupby("case", sort=False):
    g = g.sort_values("minute"); v = g["s"].values.astype(float)
    idx = pd.date_range("2020-01-01 00:00", periods=len(v), freq="1min")
    ts = pd.Series(v, index=idx); ts.index.freq = "1min"
    hand = -100 + 200*np.mean(v[:-1440] == v[1440:])
    rows.append((case, len(v), round(float(m.sri(ts)), 4), round(float(hand), 4)))
out = pd.DataFrame(rows, columns=["case", "minutes", "SRI_pyact", "SRI_formula"]); out["diff"] = (out["SRI_pyact"] - out["SRI_formula"]).round(4)
out.to_csv(out_csv, index=False); print(out.to_string(index=False))
