#!/usr/bin/env python3
"""Step 35: which QC filter changes the consensus counts, compared with removing the
same number of cells at random (pure power loss)."""
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
os.environ["N_PERM"] = "0"
src = open(os.path.join(os.path.dirname(__file__), "05_consensus.py")).read()
exec(compile(src[:src.index("t0 = time.time()")], "05", "exec"))
OUT = os.path.join(RES, "qc_raw"); rng = np.random.default_rng(8)
Q = {}
for e in EXPS:
    q = pd.read_csv(os.path.join(OUT, f"{e}_cell_qc.tsv.gz"), sep="\t").set_index("barcode")
    cells = ds[e]["cells"].cell.str.replace(".", "-", regex=False)
    Q[e] = q.reindex(cells.values).reset_index()
filters = {"low UMI (MAD)": lambda q: q.low_umi, "low genes (MAD)": lambda q: q.low_genes,
           "mito > 5%": lambda q: q.pct_mito > 5, "plastid > 5%": lambda q: q.pct_plastid > 5,
           "doublet": lambda q: q.doublet, "not in CellRanger run": lambda q: q.umi.isna()}
def run(mask):
    F = {e: fit_exp(X[e][mask[e]], Z[e][mask[e]][:, Z[e][mask[e]].std(0) > 0], Y[e][mask[e]]) for e in EXPS}
    out = {}
    for kind in ("main", "inter"):
        mt = meta(F["Exp1"][kind], F["Exp2"][kind]); out[kind] = int(consensus(F["Exp1"][kind], F["Exp2"][kind], mt).sum())
    return out
rows = []
for name, f in filters.items():
    rem = {e: f(Q[e]).fillna(False).astype(bool).values for e in EXPS}
    if name == "not in CellRanger run":
        rem = {e: Q[e].umi.isna().values for e in EXPS}
    nrem = {e: int(rem[e].sum()) for e in EXPS}
    if sum(nrem.values()) == 0:
        continue
    r = run({e: ~rem[e] for e in EXPS})
    rr = []
    for _ in range(3):
        rm = {e: np.zeros(len(rem[e]), bool) for e in EXPS}
        for e in EXPS:
            rm[e][rng.choice(len(rem[e]), nrem[e], replace=False)] = True
        rr.append(run({e: ~rm[e] for e in EXPS}))
    rows.append(dict(filter=name, removed_Exp1=nrem["Exp1"], removed_Exp2=nrem["Exp2"], main=r["main"], inter=r["inter"],
                     main_random=np.mean([x["main"] for x in rr]), inter_random=np.mean([x["inter"] for x in rr])))
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(os.path.join(OUT, "qc_filter_dissection.tsv"), sep="\t", index=False)
