#!/usr/bin/env python3
"""Step 36: permutation null (20x) for the consensus counts under strict cell QC."""
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
os.environ["N_PERM"] = "0"
src = open(os.path.join(os.path.dirname(__file__), "05_consensus.py")).read()
exec(compile(src[:src.index("t0 = time.time()")], "05", "exec"))
OUT = os.path.join(RES, "qc_raw"); rng = np.random.default_rng(21)
keep = {}
for e in EXPS:
    q = pd.read_csv(os.path.join(OUT, f"{e}_cell_qc.tsv.gz"), sep="\t").set_index("barcode")
    qq = q.reindex(ds[e]["cells"].cell.str.replace(".", "-", regex=False).values)
    bad = (qq.low_umi | qq.low_genes | (qq.pct_mito > 5) | (qq.pct_plastid > 5) | qq.doublet).fillna(True).astype(bool).values
    keep[e] = ~bad
def sub(e, perm=False):
    k = keep[e]; Xe, Ze, Ye = X[e][k], Z[e][k], Y[e][k]; Ze = Ze[:, Ze.std(0) > 0]
    if perm:
        tot = np.log1p(ds[e]["T"].sum(1).values[k]); lib = ds[e]["cells"].library.values[k]
        key = pd.Series(lib).astype(str) + "_" + pd.qcut(tot + rng.uniform(0, 1e-6, tot.size), 10, labels=False).astype(str)
        p = np.arange(len(tot))
        for kk in key.unique():
            i = np.where(key == kk)[0]; p[i] = rng.permutation(i)
        Xe = Xe[p]
    return fit_exp(Xe, Ze, Ye)
def counts(F):
    return {k: int(consensus(F["Exp1"][k], F["Exp2"][k], meta(F["Exp1"][k], F["Exp2"][k])).sum()) for k in ("main", "inter")}
real = counts({e: sub(e) for e in EXPS}); print("real", real, flush=True)
null = [counts({e: sub(e, True) for e in EXPS}) for _ in range(20)]
nd = pd.DataFrame(null); print("null mean", nd.mean().to_dict(), "max", nd.max().to_dict())
res = dict(real=real, null_mean=nd.mean().to_dict(), null_max=nd.max().to_dict(),
           emp_FDR={k: float(nd[k].mean() / max(real[k], 1)) for k in real}, cells={e: int(keep[e].sum()) for e in EXPS})
json.dump(res, open(os.path.join(OUT, "strict_qc_null.json"), "w"), indent=1); print(res)
