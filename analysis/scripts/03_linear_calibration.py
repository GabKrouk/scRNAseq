#!/usr/bin/env python3
"""Step 3: calibrate how TF presence is encoded and which covariates are
adjusted, using reproducibility between the two independent experiments as
the objective (plus a stratified permutation null for each configuration).

Model per gene g and experiment e:  y_g = a + X b_g + Z c_g + error
(log1p CP10k; HC1-robust SEs).  For every configuration we compare the
TF x gene t-statistics of Exp1 and Exp2.
"""
import os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import *

N_PERM = int(os.environ.get("N_PERM", 5))
rng = np.random.default_rng(20261004)

ds = {e: load(e) for e in EXPS}
genes, tfs = analysis_sets(ds)
print(f"{len(genes)} genes, {len(tfs)} TFs: {tfs}")
Y = {e: response(ds[e], genes) for e in EXPS}

COV = {
    "none": lambda d: covariates(d, tfs)[[c for c in covariates(d, tfs).columns if c.startswith("lib")]],
    "depth": lambda d: covariates(d, tfs)[[c for c in covariates(d, tfs).columns if c in ("log_umi", "log_ngenes") or c.startswith("lib")]],
    "depth+load": lambda d: covariates(d, tfs),
}
SCHEMES = ["ge1", "ge2", "ge3", "ge5", "log", "share"]


def metrics(f1, f2):
    t1, t2 = f1["t"], f2["t"]
    q1, q2 = bh(f1["p"]), bh(f2["p"])
    s1, s2 = q1 < 0.05, q2 < 0.05
    same = np.sign(t1) == np.sign(t2)
    rep = s1 & s2 & same
    r12 = np.mean((f2["p"][s1] < 0.05) & same[s1]) if s1.any() else np.nan
    r21 = np.mean((f1["p"][s2] < 0.05) & same[s2]) if s2.any() else np.nan
    return dict(r_t=np.corrcoef(t1.ravel(), t2.ravel())[0, 1],
                n_sig1=int(s1.sum()), n_sig2=int(s2.sum()), n_rep=int(rep.sum()),
                n_opposite=int((s1 & s2 & ~same).sum()),
                rep_rate_1to2=r12, rep_rate_2to1=r21)


rows, per_tf = [], []
for cov_name, cov_fn in COV.items():
    Z = {e: cov_fn(ds[e]).values for e in EXPS}
    for sc in SCHEMES:
        t0 = time.time()
        X = {e: encode(ds[e]["T"], tfs, sc) for e in EXPS}
        fit = {e: ols(X[e], Z[e], Y[e]) for e in EXPS}
        m = metrics(fit["Exp1"], fit["Exp2"])
        null = []
        for _ in range(N_PERM):
            fp = {}
            for e in EXPS:
                perm = strata_permutation(ds[e], rng)
                fp[e] = ols(X[e][perm], Z[e], Y[e])
            null.append(metrics(fp["Exp1"], fp["Exp2"]))
        nd = pd.DataFrame(null)
        m.update(cov=cov_name, scheme=sc,
                 null_n_rep_mean=nd.n_rep.mean(), null_n_rep_max=nd.n_rep.max(),
                 null_r_t_mean=nd.r_t.mean(), null_n_sig1=nd.n_sig1.mean(), null_n_sig2=nd.n_sig2.mean())
        m["emp_FDR_rep"] = m["null_n_rep_mean"] / max(m["n_rep"], 1)
        rows.append(m)
        for j, tf in enumerate(tfs):
            per_tf.append(dict(cov=cov_name, scheme=sc, tf=tf,
                               r_t=np.corrcoef(fit["Exp1"]["t"][j], fit["Exp2"]["t"][j])[0, 1],
                               n_rep=int(((bh(fit["Exp1"]["p"])[j] < .05) & (bh(fit["Exp2"]["p"])[j] < .05) &
                                          (np.sign(fit["Exp1"]["t"][j]) == np.sign(fit["Exp2"]["t"][j]))).sum())))
        print(f"{cov_name:11s} {sc:6s} r_t={m['r_t']:.3f} (null {m['null_r_t_mean']:.3f}) "
              f"sig1={m['n_sig1']} sig2={m['n_sig2']} rep={m['n_rep']} (null {m['null_n_rep_mean']:.1f}) "
              f"opp={m['n_opposite']} rr12={m['rep_rate_1to2']:.2f} rr21={m['rep_rate_2to1']:.2f} [{time.time()-t0:.0f}s]",
              flush=True)

pd.DataFrame(rows).to_csv(os.path.join(RES, "calibration_linear.tsv"), sep="\t", index=False)
pd.DataFrame(per_tf).to_csv(os.path.join(RES, "calibration_linear_per_tf.tsv"), sep="\t", index=False)
