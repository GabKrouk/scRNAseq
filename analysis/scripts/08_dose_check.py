#!/usr/bin/env python3
"""Step 8: dose-response check per TF.  Signature of a TF learned in one
experiment (its q<0.05 up- and down-regulated targets there) is scored in the
other experiment; the score (up minus down, residualised on covariates) is
then compared across bins of that TF's junction-read count (0, 1, 2, 3-5, 6+)
holding the other TFs fixed.  A real TF effect should show a monotone dose
response in the held-out experiment."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import *

ds = {e: load(e) for e in EXPS}
genes, tfs = analysis_sets(ds)
st = np.load(os.path.join(RES, "consensus_ge1", "full_stats.npz"))
CT = {e: pd.read_csv(os.path.join(RES, f"{e}_celltype_scores.tsv"), sep="\t").drop(columns=["cell_type", "cell"]) for e in EXPS}
Y = {e: response(ds[e], genes) for e in EXPS}
Z = {e: np.column_stack([np.ones(len(Y[e])), covariates(ds[e], tfs, flexible_load=True).values, CT[e].values]) for e in EXPS}
Yz = {e: (Y[e] - Y[e].mean(0)) / (Y[e].std(0) + 1e-9) for e in EXPS}
bins = [(0, 0), (1, 1), (2, 2), (3, 5), (6, 10 ** 6)]
rows = []
for disc, val in (("Exp1", "Exp2"), ("Exp2", "Exp1"), ("Exp1", "Exp1"), ("Exp2", "Exp2")):
    q = bh(st[f"{disc}_main_p"]); b = st[f"{disc}_main_beta"]
    T = ds[val]["T"][tfs].values
    for j, tf in enumerate(tfs):
        up = (q[j] < 0.05) & (b[j] > 0); dn = (q[j] < 0.05) & (b[j] < 0)
        if up.sum() + dn.sum() < 10:
            continue
        s = (Yz[val][:, up].mean(1) if up.any() else 0) - (Yz[val][:, dn].mean(1) if dn.any() else 0)
        # residualise on covariates and the OTHER TFs' presence
        others = (np.delete(T, j, 1) >= 1).astype(float)
        D = np.column_stack([Z[val], others])
        r = s - D @ np.linalg.lstsq(D, s, rcond=None)[0]
        row = dict(signature_from=disc, scored_in=val, tf=tf, n_up=int(up.sum()), n_down=int(dn.sum()))
        base = r[T[:, j] == 0].mean()
        for lo, hi in bins:
            m = (T[:, j] >= lo) & (T[:, j] <= hi)
            row[f"reads_{lo}-{hi if hi < 1e6 else 'max'}"] = (r[m].mean() - base) if m.sum() >= 20 else np.nan
            row[f"n_{lo}-{hi if hi < 1e6 else 'max'}"] = int(m.sum())
        # slope of score on log1p(count) among cells with >=1 read
        pos = T[:, j] >= 1
        row["dose_slope_pos_cells"] = np.polyfit(np.log1p(T[pos, j]), r[pos], 1)[0] if pos.sum() > 30 else np.nan
        rows.append(row)
df = pd.DataFrame(rows)
df.to_csv(os.path.join(RES, "dose_response_signatures.tsv"), sep="\t", index=False)
pd.set_option("display.width", 250)
cols = ["signature_from", "scored_in", "tf", "n_up", "n_down"] + [c for c in df.columns if c.startswith("reads_")] + ["dose_slope_pos_cells"]
print(df[cols].round(3).to_string(index=False))
