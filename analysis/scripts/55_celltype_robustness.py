#!/usr/bin/env python3
"""Step 55: are reproducible TF effects driven by cell identity? Refit the final models (log dose)
without the 9 cell-type marker scores and compare: main effects (standard baseline) and non-additive
interactions (load-modulated baseline)."""
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
os.environ["N_PERM"] = "0"; os.environ["SCHEME"] = "log"
src = open(os.path.join(os.path.dirname(__file__), "05_consensus.py")).read()
exec(compile(src[:src.index("t0 = time.time()")], "05", "exec"))
nct = CT["Exp1"].shape[1]
L = {e: covariates(ds[e], tfs)["log_tf_total"].values for e in EXPS}; L = {e: (v - v.mean()) / v.std() for e, v in L.items()}
def calls(F, k, q=0.05, p=0.05):
    r1, r2 = F["Exp1"][k], F["Exp2"][k]; mt = meta(r1, r2)
    return (mt["q"] < q) & (r1["p"] < p) & (r2["p"] < p) & (np.sign(r1["beta"]) == np.sign(r2["beta"])), mt
out = {}
for lab, drop in (("with_celltype", False), ("without_celltype", True)):
    Zs = {e: Z[e][:, :-nct] if drop else Z[e] for e in EXPS}
    Fm = {e: fit_exp(X[e], Zs[e], Y[e]) for e in EXPS}
    Fi = {e: fit_exp(X[e], np.column_stack([Zs[e], X[e] * L[e][:, None]]), Y[e]) for e in EXPS}
    cm, mm = calls(Fm, "main"); ci, mi = calls(Fi, "inter", 0.10, 0.10)
    out[lab] = dict(cm=cm, ci=ci, bm=mm["beta"], bi=mi["beta"])
a, b = out["with_celltype"], out["without_celltype"]
res = dict(main_with=int(a["cm"].sum()), main_without=int(b["cm"].sum()), main_overlap=int((a["cm"] & b["cm"]).sum()),
           main_beta_r=float(np.corrcoef(a["bm"].ravel(), b["bm"].ravel())[0, 1]),
           inter_with=int(a["ci"].sum()), inter_without=int(b["ci"].sum()), inter_overlap=int((a["ci"] & b["ci"]).sum()),
           inter_beta_r=float(np.corrcoef(a["bi"].ravel(), b["bi"].ravel())[0, 1]))
print(json.dumps(res, indent=1)); json.dump(res, open(os.path.join(RES, "paper", "celltype_robustness.json"), "w"), indent=1)
