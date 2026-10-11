#!/usr/bin/env python3
"""Step 44: held-out replication of dose-based interactions (log_loadmod): discover in
one experiment only (BH q<0.05), then check the other experiment (same sign, p<0.05),
compared with a permutation null of the discovery experiment's TF labels."""
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
os.environ.setdefault("N_PERM", "0")
src = open(os.path.join(os.path.dirname(__file__), "05_consensus.py")).read()
exec(compile(src[:src.index("t0 = time.time()")], "05", "exec"))
rng = np.random.default_rng(44)
Xd = {e: encode(ds[e]["T"], tfs, "log") for e in EXPS}
L = {e: covariates(ds[e], tfs)["log_tf_total"].values for e in EXPS}; L = {e: (v - v.mean()) / v.std() for e, v in L.items()}
def fit(e, perm=False):
    Xe = Xd[e][strata_permutation(ds[e], rng)] if perm else Xd[e]
    return fit_exp(Xe, np.column_stack([Z[e], Xe * L[e][:, None]]), Y[e])["inter"]
F = {e: fit(e) for e in EXPS}
out = {}
for d, v in (("Exp1", "Exp2"), ("Exp2", "Exp1")):
    disc = bh(F[d]["p"]) < 0.05
    same = np.sign(F[d]["beta"]) == np.sign(F[v]["beta"])
    rep = disc & same & (F[v]["p"] < 0.05)
    # null: permute the discovery experiment, keep validation real
    nulls = []
    for _ in range(10):
        Fp = fit(d, True); dp = bh(Fp["p"]) < 0.05
        nulls.append((int(dp.sum()), float(np.mean(np.sign(Fp["beta"][dp]) == np.sign(F[v]["beta"][dp]))) if dp.any() else np.nan))
    out[f"{d}->{v}"] = dict(discovered=int(disc.sum()), same_sign=float(same[disc].mean()) if disc.any() else None,
                            replicated_p05=int(rep.sum()), replication_rate=float(rep.sum() / max(disc.sum(), 1)),
                            null_discovered_mean=float(np.mean([n[0] for n in nulls])))
    print(d, "->", v, out[f"{d}->{v}"], flush=True)
json.dump(out, open(os.path.join(RES, "combinatorial_push", "heldout_replication.json"), "w"), indent=1)
