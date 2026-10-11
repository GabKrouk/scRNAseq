#!/usr/bin/env python3
"""Step 67: test each connection of the paper's network (Fig 5B: 461 interactions, TFa + TFb -> gate -> target)
in the Shahan et al. (2022) WT root atlas, where the TFs are ENDOGENOUS (no over-expression).

For every TF pair, cells fall in four states by endogenous detection (>= 1 raw UMI at the TF locus):
neither / TFa only / TFb only / both.  For each network target of the pair:
  y = log1p(CP10k) ~ [TFa only] + [TFb only] + [both] + stratum + covariates     (HC1 robust covariance)
  stratum = sample x cell type x developmental time group; covariates = log depth, (log depth)^2,
  log genes detected, % mitochondrial reads (atlas_common.py).
dA, dB, dAB = effects of the three states vs 'neither' (same quantities as Fig 5C-E);
non-additivity = dAB - (dA + dB) = the gate, tested per connection (BH over the 461 connections).
Gate shape in the atlas, same definitions as the paper (singles significant at p < 0.05):
  emergent    neither TF alone acts, the pair does
  synergy     the pair goes beyond the sum, in the direction of the single-TF effect(s)
  competition the pair falls back towards 'neither' (deviation opposite to the single-TF effects)
  reversal    the pair goes opposite to the single-TF effect(s)
"""
import os, sys, json
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from common import RES, bh
from atlas_common import load_atlas, AGI

OUT = os.path.join(RES, "atlas"); os.makedirs(OUT, exist_ok=True)
NET = pd.read_csv(os.path.join(RES, "paper", "fig5d_candidates.tsv"), sep="\t")
NET["gate"] = NET["shape"].replace({"cancellation": "competition"})
NET["a"] = NET.pair.str.split("-").str[0]; NET["b"] = NET.pair.str.split("-").str[1]
A = load_atlas()
n, tilde, tot = A["n"], A["tilde"], A["tot"]
need = sorted(set(NET.target) | {AGI[t] for t in set(NET.a) | set(NET.b)})
missing = [g for g in need if g not in A["genes"]]
print("not in atlas:", missing)
need = [g for g in need if g in A["genes"]]
R = A["M"][A["genes"].get_indexer(need), :].tocsr()
row = {g: np.asarray(R[i].todense()).ravel().astype(np.float64) for i, g in enumerate(need)}
del R
det = {t: row[AGI[t]] >= 1 for t in set(NET.a) | set(NET.b)}

rows, plot_rows = [], []
for p, sub in NET.groupby("pair"):
    a, b = p.split("-")
    sa, sb = det[a], det[b]
    S3 = np.column_stack([sa & ~sb, ~sa & sb, sa & sb]).astype(np.float64)
    ncell = dict(n_neither=int((~sa & ~sb).sum()), n_A=int(S3[:, 0].sum()), n_B=int(S3[:, 1].sum()), n_both=int(S3[:, 2].sum()))
    X = tilde(S3)
    XtXi = np.linalg.inv(X.T @ X)
    for r in sub.itertuples():
        if r.target not in row:
            continue
        y = tilde(np.log1p(row[r.target] / tot * 1e4)[:, None])[:, 0]
        beta = XtXi @ (X.T @ y)
        e = y - X @ beta
        V = XtXi @ ((X * (e ** 2)[:, None]).T @ X) @ XtXi * n / (A["DF"] - 3)
        se = np.sqrt(np.diag(V))
        c = np.array([-1.0, -1.0, 1.0]); na = c @ beta; na_se = np.sqrt(c @ V @ c)
        out = dict(pair=p, a=a, b=b, target=r.target, gate_paper=r.gate, tier=r.tier,
                   paper_dA=r.dA, paper_dB=r.dB, paper_dAB=r.dAB, paper_nonadd=r.nonadd, **ncell,
                   dA=beta[0], dB=beta[1], dAB=beta[2], se_dA=se[0], se_dB=se[1], se_dAB=se[2],
                   p_dA=2 * stats.norm.sf(abs(beta[0] / se[0])), p_dB=2 * stats.norm.sf(abs(beta[1] / se[1])),
                   p_dAB=2 * stats.norm.sf(abs(beta[2] / se[2])),
                   nonadd=na, se_nonadd=na_se, z_nonadd=na / na_se, p_nonadd=2 * stats.norm.sf(abs(na / na_se)),
                   target_det_frac=float((row[r.target] > 0).mean()))
        rows.append(out)
        for st, est, s_ in (("neither", 0.0, 0.0), ("A", beta[0], se[0]), ("B", beta[1], se[1]), ("both", beta[2], se[2])):
            plot_rows.append(dict(pair=p, target=r.target, state=st, est=est, se=s_))
    print(p, ncell, flush=True)

T = pd.DataFrame(rows)
T["q_nonadd"] = bh(T.p_nonadd.values)
T["same_sign"] = np.sign(T.nonadd) == np.sign(T.paper_nonadd)
T["status"] = np.where(T.q_nonadd >= 0.05, "not significant", np.where(T.same_sign, "recovered", "opposite"))
# atlas gate shape (only meaningful when the gate is significant)
sA, sB = T.p_dA < 0.05, T.p_dB < 0.05
s_single = np.sign(np.where(sA, T.dA, 0) + np.where(sB, T.dB, 0))
T["gate_atlas"] = np.select(
    [T.q_nonadd >= 0.05, ~sA & ~sB, np.sign(T.dAB) == -s_single, np.sign(T.nonadd) == s_single],
    ["none", "emergent", "reversal", "synergy"], default="competition")
T["same_gate"] = T.gate_atlas == T.gate_paper
T.to_csv(os.path.join(OUT, "network_connections_atlas.tsv"), sep="\t", index=False)
pd.DataFrame(plot_rows).to_csv(os.path.join(OUT, "network_connections_states.tsv"), sep="\t", index=False)

summ = dict(
    n_connections=len(NET), n_tested=len(T), n_pairs=int(T.pair.nunique()), n_targets=int(T.target.nunique()),
    status=T.status.value_counts().to_dict(),
    same_sign_all=float(T.same_sign.mean()),
    binom_p_same_sign_among_significant=float(stats.binomtest(int((T.status == "recovered").sum()),
                                                              int((T.status != "not significant").sum()), 0.5).pvalue),
    recovered_same_gate=int(((T.status == "recovered") & T.same_gate).sum()),
    by_gate=T.groupby("gate_paper").status.value_counts().unstack(fill_value=0).to_dict("index"),
    single_TF_effects=dict(
        same_sign_dA=float((np.sign(T.dA) == np.sign(T.paper_dA)).mean()),
        same_sign_dB=float((np.sign(T.dB) == np.sign(T.paper_dB)).mean()),
        rho_dA=float(stats.spearmanr(T.dA, T.paper_dA).correlation),
        rho_dB=float(stats.spearmanr(T.dB, T.paper_dB).correlation),
        rho_dAB=float(stats.spearmanr(T.dAB, T.paper_dAB).correlation)),
    rho_nonadd=float(stats.spearmanr(T.nonadd, T.paper_nonadd).correlation))
json.dump(summ, open(os.path.join(OUT, "network_connections_summary.json"), "w"), indent=1, default=float)
print(json.dumps(summ, indent=1, default=float))
