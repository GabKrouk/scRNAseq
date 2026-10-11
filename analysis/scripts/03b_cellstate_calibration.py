#!/usr/bin/env python3
"""Step 3b: does adjusting for cell identity / cell state improve cross-
experiment reproducibility?  Cell state = top-k principal components of the
endogenous transcriptome (computed per experiment, after removing depth/load
covariates).  We also check how much of each PC is explained by the TF design
(if PCs carried TF effects, adjusting for them would remove true signal).
"""
import os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import *

N_PERM = int(os.environ.get("N_PERM", 5))
rng = np.random.default_rng(7)
ds = {e: load(e) for e in EXPS}
genes, tfs = analysis_sets(ds)
Y = {e: response(ds[e], genes) for e in EXPS}
Z0 = {e: covariates(ds[e], tfs) for e in EXPS}

# PCs of covariate-residualised, standardised expression (per experiment)
PCS, KMAX = {}, 40
for e in EXPS:
    D = np.column_stack([np.ones(len(Y[e])), Z0[e].values])
    R = Y[e] - D @ np.linalg.lstsq(D, Y[e], rcond=None)[0]
    R = R / (R.std(0) + 1e-9)
    U, S, Vt = np.linalg.svd(R, full_matrices=False)
    PCS[e] = U[:, :KMAX] * S[:KMAX]
    ve = S ** 2 / (S ** 2).sum()
    print(e, "variance explained by PCs 1-10:", (ve[:10] * 100).round(2))
    for sc in ("ge1", "share"):
        X = encode(ds[e]["T"], tfs, sc)
        Xd = np.column_stack([np.ones(len(X)), X])
        r2 = []
        for k in range(10):
            y = PCS[e][:, k]
            res = y - Xd @ np.linalg.lstsq(Xd, y, rcond=None)[0]
            r2.append(1 - res.var() / y.var())
        print(f"  R2 of PC1-10 explained by the {sc} TF design:", np.round(r2, 3))
np.save(os.path.join(DATA, "Exp1_pcs.npy"), PCS["Exp1"]); np.save(os.path.join(DATA, "Exp2_pcs.npy"), PCS["Exp2"])


def metrics(f1, f2):
    t1, t2 = f1["t"], f2["t"]
    q1, q2 = bh(f1["p"]), bh(f2["p"])
    s1, s2 = q1 < 0.05, q2 < 0.05
    same = np.sign(t1) == np.sign(t2)
    # TF-specific component: remove per-gene mean t across TFs
    c1, c2 = t1 - t1.mean(0), t2 - t2.mean(0)
    return dict(r_t=np.corrcoef(t1.ravel(), t2.ravel())[0, 1],
                r_t_specific=np.corrcoef(c1.ravel(), c2.ravel())[0, 1],
                r_common=np.corrcoef(t1.mean(0), t2.mean(0))[0, 1],
                n_sig1=int(s1.sum()), n_sig2=int(s2.sum()),
                n_rep=int((s1 & s2 & same).sum()), n_opp=int((s1 & s2 & ~same).sum()))


rows = []
for sc in ("ge1", "ge2", "ge3", "log", "share"):
    X = {e: encode(ds[e]["T"], tfs, sc) for e in EXPS}
    for k in (0, 5, 10, 20, 40):
        Z = {e: np.column_stack([Z0[e].values, PCS[e][:, :k]]) for e in EXPS}
        f = {e: ols(X[e], Z[e], Y[e]) for e in EXPS}
        m = metrics(f["Exp1"], f["Exp2"])
        nl = []
        for _ in range(N_PERM):
            fp = {e: ols(X[e][strata_permutation(ds[e], rng)], Z[e], Y[e]) for e in EXPS}
            nl.append(metrics(fp["Exp1"], fp["Exp2"]))
        nl = pd.DataFrame(nl)
        m.update(scheme=sc, n_pcs=k, null_r_t=nl.r_t.mean(), null_n_rep=nl.n_rep.mean(),
                 null_n_sig1=nl.n_sig1.mean(), null_n_sig2=nl.n_sig2.mean())
        rows.append(m)
        print(sc, k, {a: (round(b, 3) if isinstance(b, float) else b) for a, b in m.items()}, flush=True)
pd.DataFrame(rows).to_csv(os.path.join(RES, "calibration_cellstate.tsv"), sep="\t", index=False)
