#!/usr/bin/env python3
"""Step 17: Fig 6A/6B genes (qPCR-validated DoubleTARGET targets, Sup Table 3)
re-examined in Exp1, Exp2 and both combined with the four-state model of
step 13 (neither / ERF5 only / partner only / both), for ERF5-CRF4 and ERF5-NAC4."""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from common import *

GENES = ["AT1G65845", "AT1G13390", "AT1G19210", "AT1G55920", "AT2G26530", "AT2G38870", "AT2G40140",
         "AT3G29034", "AT3G46080", "AT3G52450", "AT4G21920", "AT4G33666", "AT4G36500", "AT5G07165",
         "AT5G22520", "AT5G24110", "AT5G39670", "AT5G43650", "AT5G62480"]
PAIRS = [("ERF5", "CRF4"), ("ERF5", "NAC4")]
OUT = os.path.join(RES, "fig6"); os.makedirs(OUT, exist_ok=True)
ds = {e: load(e) for e in EXPS}
_, tfs = analysis_sets(ds)
CT = {e: pd.read_csv(os.path.join(RES, f"{e}_celltype_scores.tsv"), sep="\t").drop(columns=["cell_type", "cell"]) for e in EXPS}
present = [g for g in GENES if all(g in ds[e]["genes"] for e in EXPS)]
print("genes found in both matrices:", len(present), "missing:", sorted(set(GENES) - set(present)))
det = {e: pd.Series((ds[e]["G"][pd.Index(ds[e]["genes"]).get_indexer(present)] > 0).mean(1).A1, index=present) for e in EXPS}
Y = {e: response(ds[e], present) for e in EXPS}
Z = {e: np.column_stack([covariates(ds[e], tfs, flexible_load=True).values, CT[e].values]) for e in EXPS}
X = {e: encode(ds[e]["T"], tfs, "ge1") for e in EXPS}
tfi = {t: i for i, t in enumerate(tfs)}
means, tests = [], []
for ta, tb in PAIRS:
    a, b = tfi[ta], tfi[tb]
    pooled = []
    for e in EXPS:
        xa, xb = X[e][:, a], X[e][:, b]
        S = np.column_stack([xa * (1 - xb), (1 - xa) * xb, xa * xb])
        D = np.column_stack([np.ones(len(xa)), S, np.delete(X[e], [a, b], 1), Z[e]])
        A = np.linalg.pinv(D.T @ D) @ D.T; B = A @ Y[e]; E = Y[e] - D @ B
        h = np.einsum("ij,ji->i", D, A); W = (E / (1 - h)[:, None]) ** 2
        Lc = np.array([[0, -1, -1, 1]], float)  # interaction = both - a - b (on S coefs)
        As = A[1:4]
        inter = (Lc[:, 1:] @ B[1:4])[0]
        vint = ((Lc[:, 1:] @ As) ** 2) @ W
        adj = Y[e] - D[:, 4:] @ B[4:]
        state = np.select([xa * xb == 1, xa == 1, xb == 1], ["both", "a", "b"], "neither")
        for j, g in enumerate(present):
            ref = adj[state == "neither", j].mean()
            for s in ("neither", "a", "b", "both"):
                v = adj[state == s, j] - ref
                means.append(dict(pair=f"{ta}-{tb}", TFa=ta, TFb=tb, target=g, exp=e, state=s, mean=v.mean(),
                                  ci=1.96 * v.std() / np.sqrt(len(v)), n=len(v)))
            pooled.append(pd.DataFrame({"target": g, "state": state, "delta": adj[:, j] - ref}))
            tests.append(dict(pair=f"{ta}-{tb}", target=g, exp=e, interaction=inter[j],
                              p_interaction=2 * stats.norm.sf(abs(inter[j]) / np.sqrt(vint[0, j])),
                              effect_a=B[1, j], effect_b=B[2, j], effect_both=B[3, j]))
    P = pd.concat(pooled)
    for (g, s), v in P.groupby(["target", "state"]).delta:
        means.append(dict(pair=f"{ta}-{tb}", TFa=ta, TFb=tb, target=g, exp="Combined", state=s, mean=v.mean(),
                          ci=1.96 * v.std() / np.sqrt(len(v)), n=len(v)))
M = pd.DataFrame(means); T = pd.DataFrame(tests)
# meta interaction p
w = T.pivot_table(index=["pair", "target"], columns="exp", values=["interaction", "p_interaction"])
for e in EXPS:
    T.loc[T.exp == e, "se"] = (T.interaction / stats.norm.isf(T.p_interaction / 2)).abs()
m = T.pivot_table(index=["pair", "target"], columns="exp", values=["interaction", "se", "p_interaction"])
wt = 1 / m["se"] ** 2
m[("interaction", "meta")] = (m["interaction"] * wt).sum(1) / wt.sum(1)
m[("p_interaction", "meta")] = 2 * stats.norm.sf(abs(m[("interaction", "meta")]) * np.sqrt(wt.sum(1)))
m.columns = [f"{a}_{b}" for a, b in m.columns]
m = m.reset_index()
m["same_sign"] = np.sign(m.interaction_Exp1) == np.sign(m.interaction_Exp2)
m["detect_Exp1"] = m.target.map(det["Exp1"]); m["detect_Exp2"] = m.target.map(det["Exp2"])
M.to_csv(os.path.join(OUT, "fig6_state_means.tsv"), sep="\t", index=False)
m.to_csv(os.path.join(OUT, "fig6_interaction_tests.tsv"), sep="\t", index=False)
pd.set_option("display.width", 250)
print(m[["pair", "target", "interaction_Exp1", "p_interaction_Exp1", "interaction_Exp2", "p_interaction_Exp2",
         "interaction_meta", "p_interaction_meta", "same_sign", "detect_Exp1", "detect_Exp2"]].round(3).to_string(index=False))
g = M[(M.target == "AT1G65845") & (M.pair == "ERF5-CRF4")]
print(g.pivot(index="state", columns="exp", values="mean").round(3))

# ---- Fig 6A analogue: AT1G65845, depth-normalised only (as in the paper) vs fully adjusted
rows = []
for thrC in (1, 3):
    for e in EXPS:
        d = ds[e]; gi_ = d["genes"].index("AT1G65845")
        y = d["G"][gi_, :].toarray().ravel() / d["G"].sum(0).A1 * 1e4
        a = d["T"].ERF5.values >= 1; b = d["T"].CRF4.values >= thrC
        ref = y[~a & ~b].mean()
        for s, msk in (("neither", ~a & ~b), ("a", a & ~b), ("b", ~a & b), ("both", a & b)):
            v = y[msk] / ref * 100
            rows.append(dict(CRF4_threshold=thrC, exp=e, state=s, mean=v.mean(), se=v.std() / np.sqrt(msk.sum()), n=int(msk.sum())))
pd.DataFrame(rows).to_csv(os.path.join(OUT, "fig6A_AT1G65845_depthnorm.tsv"), sep="\t", index=False)
