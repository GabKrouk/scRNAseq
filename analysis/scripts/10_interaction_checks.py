#!/usr/bin/env python3
"""Step 10: characterise the reproducible TF x TF interactions.
(a) super-additive (interaction has the same sign as the summed main effects)
    vs sub-additive (opposite sign: the pair does less than the sum of parts);
(b) dose sensitivity: for the pairs with most consensus interactions, refit
    adding each TF's own dose (log1p junction reads) and its partner's, so an
    'interaction' that only reflects one TF diluting the other's dose is
    removed (the concern Reviewer 3 raised for DoubleTARGET);
    A second refit also lets each TF's main effect scale with GR-TF load
    (x_a * load, x_b * load), so a pair 'interaction' that is really one TF
    acting more strongly in heavily transfected cells is removed;
(c) the ERF5-CRF4 and ERF5-NAC4 pairs from the manuscript."""
import os, sys, itertools, json
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from common import *

OUT = os.path.join(RES, "consensus_ge1")
st = np.load(os.path.join(OUT, "full_stats.npz"))
genes, tfs, pairs = list(st["genes"]), list(st["tfs"]), list(st["pairs"])
tfi = {t: i for i, t in enumerate(tfs)}


def meta(b1, s1, b2, s2):
    w1, w2 = 1 / s1 ** 2, 1 / s2 ** 2
    b = (b1 * w1 + b2 * w2) / (w1 + w2)
    return b, 2 * stats.norm.sf(np.abs(b) * np.sqrt(w1 + w2))


bm, pm = meta(st["Exp1_main_beta"], st["Exp1_main_se"], st["Exp2_main_beta"], st["Exp2_main_se"])
bi, pi_ = meta(st["Exp1_inter_beta"], st["Exp1_inter_se"], st["Exp2_inter_beta"], st["Exp2_inter_se"])
cons_i = (bh(pi_) < .05) & (st["Exp1_inter_p"] < .05) & (st["Exp2_inter_p"] < .05) & \
         (np.sign(st["Exp1_inter_beta"]) == np.sign(st["Exp2_inter_beta"]))
rows = []
for k, pr in enumerate(pairs):
    a, b = (tfi[x] for x in pr.split("-"))
    for g in np.where(cons_i[k])[0]:
        add = bm[a, g] + bm[b, g]
        rows.append(dict(pair=pr, target=genes[g], beta_TFa=bm[a, g], beta_TFb=bm[b, g], beta_interaction=bi[k, g],
                         q_interaction_meta=bh(pi_)[k, g],
                         type="super-additive" if np.sign(bi[k, g]) == np.sign(add) else "sub-additive"))
ti = pd.DataFrame(rows)
ti.to_csv(os.path.join(OUT, "consensus_interactions_classified.tsv"), sep="\t", index=False)
print("consensus interactions by type:", ti["type"].value_counts().to_dict())
print(ti.groupby("pair")["type"].value_counts().unstack(fill_value=0).assign(
    total=lambda d: d.sum(1)).sort_values("total", ascending=False).head(15).to_string())

# (b) dose sensitivity for top pairs
ds = {e: load(e) for e in EXPS}
CT = {e: pd.read_csv(os.path.join(RES, f"{e}_celltype_scores.tsv"), sep="\t").drop(columns=["cell_type", "cell"]) for e in EXPS}
Y = {e: response(ds[e], genes) for e in EXPS}
Zb = {e: np.column_stack([covariates(ds[e], tfs, flexible_load=True).values, CT[e].values]) for e in EXPS}
X = {e: encode(ds[e]["T"], tfs, "ge1") for e in EXPS}
top = ti["pair"].value_counts().head(10).index.tolist() + ["CRF4-ERF5", "ERF5-NAC4"]
# global BH cut-off used for the original consensus (so refits use the same bar)
q_all = bh(pi_)
p_cut = np.nanmax(np.where(q_all < 0.05, pi_, np.nan))
res = []
for pr in dict.fromkeys(top):
    a, b = (tfi[x] for x in pr.split("-"))
    k = pairs.index(pr)
    out = dict(pair=pr, consensus_original=int(cons_i[k].sum()))
    for label in ("dose", "dose+load_modulation"):
        fits = {}
        for e in EXPS:
            extra = [np.log1p(ds[e]["T"][[tfs[a], tfs[b]]].values)]
            if label == "dose+load_modulation":
                L = covariates(ds[e], tfs)["log_tf_total"].values
                L = (L - L.mean()) / L.std()
                # each TF's effect allowed to scale with GR-TF load
                extra.append(np.column_stack([X[e][:, a] * L, X[e][:, b] * L]))
            Xi = np.column_stack([X[e], X[e][:, a] * X[e][:, b]])
            fits[e] = ols(Xi, np.column_stack([Zb[e]] + extra), Y[e])
        b1, s1, p1 = fits["Exp1"]["beta"][-1], fits["Exp1"]["se"][-1], fits["Exp1"]["p"][-1]
        b2, s2, p2 = fits["Exp2"]["beta"][-1], fits["Exp2"]["se"][-1], fits["Exp2"]["p"][-1]
        bmd, pmd = meta(b1, s1, b2, s2)
        c = (pmd <= p_cut) & (p1 < .05) & (p2 < .05) & (np.sign(b1) == np.sign(b2))
        out[f"retained_{label}"] = int((cons_i[k] & c).sum())
        out[f"new_{label}"] = int((~cons_i[k] & c).sum())
        out[f"corr_beta_{label}"] = float(np.corrcoef(bi[k], bmd)[0, 1])
    res.append(out)
    print(out, flush=True)
pd.DataFrame(res).to_csv(os.path.join(OUT, "interaction_dose_sensitivity.tsv"), sep="\t", index=False)
