#!/usr/bin/env python3
"""Step 9: how recovery of reproducible TF->target edges depends on the number
of cells.  One experiment is randomly down-sampled (3 replicates per size),
the other is kept whole; we count consensus main-effect edges (same criteria
as step 5) overall and per TF, together with the number of TF-positive cells
behind each estimate."""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from common import *

rng = np.random.default_rng(5)
ds = {e: load(e) for e in EXPS}
genes, tfs = analysis_sets(ds)
CT = {e: pd.read_csv(os.path.join(RES, f"{e}_celltype_scores.tsv"), sep="\t").drop(columns=["cell_type", "cell"]) for e in EXPS}
Y = {e: response(ds[e], genes) for e in EXPS}
Z = {e: np.column_stack([covariates(ds[e], tfs, flexible_load=True).values, CT[e].values]) for e in EXPS}
X = {e: encode(ds[e]["T"], tfs, os.environ.get("SCHEME", "log")) for e in EXPS}
full = {e: ols(X[e], Z[e], Y[e]) for e in EXPS}


def consensus(a, b):
    w1, w2 = 1 / a["se"] ** 2, 1 / b["se"] ** 2
    bm = (a["beta"] * w1 + b["beta"] * w2) / (w1 + w2)
    pm = 2 * stats.norm.sf(np.abs(bm) * np.sqrt(w1 + w2))
    return (bh(pm) < 0.05) & (a["p"] < 0.05) & (b["p"] < 0.05) & (np.sign(a["beta"]) == np.sign(b["beta"]))


rows = []
for sub, keep in (("Exp2", "Exp1"), ("Exp1", "Exp2")):
    n_all = len(Y[sub])
    sizes = [s for s in (500, 1000, 1500, 1930, 3000, 4500, 5811) if s <= n_all]
    for n in sizes:
        for rep in range(3 if n < n_all else 1):
            idx = rng.choice(n_all, n, replace=False) if n < n_all else np.arange(n_all)
            Zs = Z[sub][idx]
            Zs = Zs[:, Zs.std(0) > 0]                         # drop constant library dummies
            f = ols(X[sub][idx], Zs, Y[sub][idx])
            c = consensus(f, full[keep]) if sub == "Exp2" else consensus(full[keep], f)
            row = dict(subsampled=sub, n_cells=n, rep=rep, n_consensus=int(c.sum()),
                       n_sig_subsampled=int((bh(f["p"]) < 0.05).sum()))
            for j, tf in enumerate(tfs):
                row[f"cons_{tf}"] = int(c[j].sum())
                row[f"poscells_{tf}"] = int((X[sub][idx, j] > 0).sum())
            rows.append(row)
            print(sub, n, rep, row["n_consensus"], flush=True)
df = pd.DataFrame(rows)
df.to_csv(os.path.join(RES, "power_vs_cells.tsv"), sep="\t", index=False)
