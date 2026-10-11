#!/usr/bin/env python3
"""Step 13: four-state view of every reproducible TFa x TFb -> target interaction.

For each consensus interaction (step 5) the target is refitted in each
experiment with a cell-means coding of the pair:
    y = c + e_a*[a only] + e_b*[b only] + e_ab*[a and b] + other 21 TFs + Z
so e_a, e_b, e_ab are the effects of "TFa alone", "TFb alone" and "both"
relative to "neither" (same covariates as everywhere else, HC3 SEs).
Experiments are combined by inverse-variance meta-analysis, including the
contrasts both-vs-a-alone and both-vs-b-alone (full HC3 covariance).

A target is a *proper combinatorial* case (expression of the pair cannot be
predicted from either TF alone) when:
  * 'both' differs from 'a alone' AND from 'b alone' with p<0.05 in EACH
    experiment and the same direction; 'both' and the interaction have the
    same sign in both experiments; single-TF effects do not contradict, and
  * the interaction is larger than the largest single-TF effect
    (|e_ab - e_a - e_b| > max(|e_a|, |e_b|)).
Shapes:
  AND-gate     neither single TF has an effect (p>0.05), both together do
  reversal     'both' goes opposite to the significant single-TF effect(s)
  cancellation single TF(s) act, 'both' returns to the 'neither' level
  synergy      singles act in the same direction, 'both' goes far beyond
Writes per-target tables and adjusted single-cell values for plotting.
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from common import *

OUT = os.path.join(RES, "combinatorial" + ("_nTFcat" if os.environ.get("NTF_CAT", "0") == "1" else "")); os.makedirs(OUT, exist_ok=True)
ti = pd.read_csv(os.path.join(RES, "consensus_ge1", "consensus_interactions_classified.tsv"), sep="\t")
ds = {e: load(e) for e in EXPS}
genes, tfs = analysis_sets(ds)
gi = {g: i for i, g in enumerate(genes)}; tfi = {t: i for i, t in enumerate(tfs)}
CT = {e: pd.read_csv(os.path.join(RES, f"{e}_celltype_scores.tsv"), sep="\t").drop(columns=["cell_type", "cell"]) for e in EXPS}
Y = {e: response(ds[e], genes) for e in EXPS}
Z = {e: np.column_stack([covariates(ds[e], tfs, flexible_load=True).values, CT[e].values]) for e in EXPS}
X = {e: encode(ds[e]["T"], tfs, "ge1") for e in EXPS}


NTF_CAT = os.environ.get("NTF_CAT", "0") == "1"


def fit_states(e, a, b, gidx):
    Xe = X[e]; xa, xb = Xe[:, a], Xe[:, b]
    S = np.column_stack([xa * (1 - xb), (1 - xa) * xb, xa * xb])          # a only, b only, both
    others = np.delete(Xe, [a, b], axis=1)
    n = len(xa)
    Ze = Z[e]
    if NTF_CAT:   # robustness: one dummy per number of OTHER TFs detected (0..21)
        k = others.sum(1).astype(int)
        Ze = np.column_stack([Ze] + [(k == v).astype(float) for v in np.unique(k)[1:]])
    D = np.column_stack([np.ones(n), S, others, Ze])
    Yg = Y[e][:, gidx]
    XtX = np.linalg.pinv(D.T @ D); A = XtX @ D.T
    B = A @ Yg; E = Yg - D @ B
    h = np.einsum("ij,ji->i", D, A)
    W = (E / (1 - h)[:, None]) ** 2
    As = A[1:4]
    cov = np.empty((3, 3, len(gidx)))
    for i in range(3):
        for j in range(i, 3):
            cov[i, j] = cov[j, i] = (As[i] * As[j]) @ W
    # covariate- and other-TF-adjusted single-cell values (state effects kept)
    adj = Yg - D[:, 4:] @ B[4:] + B[0]       # remove others + Z, keep intercept + states
    state = np.select([xa * xb == 1, xa == 1, xb == 1], ["both", "a", "b"], "neither")
    return B[1:4], cov, adj, state


def contrasts(beta, cov):
    """effects: a, b, both, both-a, both-b, interaction (both-a-b), with variances."""
    L = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 1], [0, -1, 1], [-1, -1, 1]], float)
    est = L @ beta
    var = np.einsum("ki,ijg,kj->kg", L, cov, L)
    return est, var


names = ["a", "b", "both", "both_vs_a", "both_vs_b", "interaction"]
rows, cells = [], []
for pr, sub in ti.groupby("pair"):
    ta, tb = pr.split("-"); a, b = tfi[ta], tfi[tb]
    gidx = [gi[g] for g in sub.target]
    res = {}
    for e in EXPS:
        beta, cov, adj, state = fit_states(e, a, b, gidx)
        res[e] = contrasts(beta, cov)
        for k, g in enumerate(sub.target):
            cells.append(pd.DataFrame({"pair": pr, "target": g, "exp": e, "state": state,
                                       "adj_expr": adj[:, k].astype(np.float32)}))
        res[e + "_n"] = {s: int((state == s).sum()) for s in ("neither", "a", "b", "both")}
    for k, g in enumerate(sub.target):
        r = dict(pair=pr, TFa=ta, TFb=tb, target=g)
        for e in EXPS:
            for s, v in res[e + "_n"].items():
                r[f"n_{s}_{e}"] = v
        for j, nm in enumerate(names):
            e1, v1 = res["Exp1"][0][j, k], res["Exp1"][1][j, k]
            e2, v2 = res["Exp2"][0][j, k], res["Exp2"][1][j, k]
            w1, w2 = 1 / v1, 1 / v2
            em = (e1 * w1 + e2 * w2) / (w1 + w2); sm = np.sqrt(1 / (w1 + w2))
            r.update({f"{nm}_Exp1": e1, f"p_{nm}_Exp1": 2 * stats.norm.sf(abs(e1) / np.sqrt(v1)),
                      f"{nm}_Exp2": e2, f"p_{nm}_Exp2": 2 * stats.norm.sf(abs(e2) / np.sqrt(v2)),
                      f"{nm}_meta": em, f"p_{nm}_meta": 2 * stats.norm.sf(abs(em) / sm)})
        rows.append(r)

df = pd.DataFrame(rows)
sig = lambda nm: df[f"p_{nm}_meta"] < 0.05
sig_each = lambda nm: (df[f"p_{nm}_Exp1"] < 0.05) & (df[f"p_{nm}_Exp2"] < 0.05)
same = lambda nm: np.sign(df[f"{nm}_Exp1"]) == np.sign(df[f"{nm}_Exp2"])
ea, eb, eab, inter = df["a_meta"], df["b_meta"], df["both_meta"], df["interaction_meta"]
# 'both' must differ from each single TF in EACH experiment, same direction
differs_both = sig_each("both_vs_a") & sig_each("both_vs_b") & same("both_vs_a") & same("both_vs_b")
dominant = inter.abs() > np.maximum(ea.abs(), eb.abs())
# single-TF effects must not contradict each other between experiments
singles_ok = ~((df["p_a_Exp1"] < .05) & (df["p_a_Exp2"] < .05) & ~same("a")) & \
             ~((df["p_b_Exp1"] < .05) & (df["p_b_Exp2"] < .05) & ~same("b"))
df["proper_combinatorial"] = differs_both & dominant & same("interaction") & same("both") & singles_ok
single_sig_a, single_sig_b = sig("a"), sig("b")
sgn_singles = np.sign(np.where(single_sig_a, ea, 0) + np.where(single_sig_b, eb, 0))
shape = np.where(~single_sig_a & ~single_sig_b, "AND-gate",
         np.where(np.sign(eab) == -sgn_singles, np.where(sig("both"), "reversal", "cancellation"),
         np.where(~sig("both"), "cancellation", "synergy")))
df["shape"] = np.where(df["proper_combinatorial"], shape, "additive-like / one TF dominant")
df["score"] = -np.log10(df["p_interaction_meta"].clip(1e-300)) * inter.abs()
df = df.sort_values(["proper_combinatorial", "score"], ascending=[False, False])
df.to_csv(os.path.join(OUT, "four_state_effects_all_consensus.tsv"), sep="\t", index=False)
prop = df[df.proper_combinatorial]
prop.to_csv(os.path.join(OUT, "proper_combinatorial_targets.tsv"), sep="\t", index=False)
print("consensus interactions:", len(df), "| unique targets:", df.target.nunique())
print("proper combinatorial:", len(prop), "| unique targets:", prop.target.nunique(), "| pairs:", prop.pair.nunique())
print(prop["shape"].value_counts().to_string())
print(prop.groupby("pair").size().sort_values(ascending=False).head(15).to_string())
cells = pd.concat(cells)
keep = set(zip(df.pair, df.target))
cells.to_csv(os.path.join(OUT, "cells_adjusted_expression.tsv.gz"), sep="\t", index=False)
