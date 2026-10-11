#!/usr/bin/env python3
"""Step 51: bulk ERF5 x CRF4 DoubleTARGET (Fig 6, ex-Fig 7) - interaction term with multiple-testing
control (R3 #7). Model per gene on log2 CPM: ~ experiment + ERF5 + CRF4 + ERF5:CRF4, with limma-style
empirical-Bayes variance moderation; BH q-values; and an empirical null from permuting construct labels
within each experiment (how many genes reach p<0.05 when no true interaction exists)."""
import os, sys, json
import numpy as np, pandas as pd
from scipy import stats, special
sys.path.insert(0, os.path.dirname(__file__))
from common import RES, bh
W = "/Users/gabkrouk/Library/CloudStorage/Dropbox/Science/Papers/Emanuel New TARGET paper/Figure Will"
OUT = os.path.join(RES, "bulk_bridge")
C = pd.read_csv(os.path.join(W, "counts.txt"), sep=" ", index_col=0).drop(columns=["E2_EC.2", "E2_E.1"])
cpm = C / C.sum() * 1e6
L = np.log2(cpm[(cpm > 1).sum(1) >= 3] + 1)
s = pd.Series(L.columns); exp2 = s.str.startswith("E2").values
con = s.str.extract(r"_(\w+)\.")[0].values
Y = L.values.T
rng = np.random.default_rng(51)


def fit(con):
    E = np.isin(con, ["E", "EC"]).astype(float); Cc = np.isin(con, ["C", "EC"]).astype(float)
    D = np.column_stack([np.ones(len(con)), exp2, E, Cc, E * Cc])
    B = np.linalg.lstsq(D, Y, rcond=None)[0]; R = Y - D @ B
    d = len(con) - D.shape[1]; s2 = (R ** 2).sum(0) / d
    # limma eBayes: fit scaled inverse-chi2 prior to s2 (method of moments on log s2)
    z = np.log(s2[s2 > 0]); e = z - special.digamma(d / 2) + np.log(d / 2)
    v = np.var(e, ddof=1) - special.polygamma(1, d / 2)
    d0 = 2 * _inv_trigamma(v) if v > 0 else np.inf
    s02 = np.exp(np.mean(e) + special.digamma(d0 / 2) - np.log(d0 / 2)) if np.isfinite(d0) else np.exp(np.mean(e))
    s2m = (d0 * s02 + d * s2) / (d0 + d) if np.isfinite(d0) else np.full_like(s2, s02)
    cv = np.linalg.inv(D.T @ D)[4, 4]
    t = B[4] / np.sqrt(s2m * cv); dft = d + (d0 if np.isfinite(d0) else 1e6)
    return B[4], 2 * stats.t.sf(np.abs(t), dft), 2 * stats.t.sf(np.abs(B[4] / np.sqrt(s2 * cv)), d)


def _inv_trigamma(x):
    y = 0.5 + 1 / x
    for _ in range(50):
        tri = special.polygamma(1, y); y = y + tri * (1 - tri / x) / special.polygamma(2, y)
    return y


b, p_mod, p_ols = fit(con)
q = bh(p_mod)
null = []
for i in range(200):
    pc = con.copy()
    for ex in (False, True):
        idx = np.where(exp2 == ex)[0]; pc[idx] = rng.permutation(pc[idx])
    null.append(int((fit(pc)[1] < 0.05).sum()))
will = None
try:
    import xlsx_min
except Exception:
    pass
res = dict(n_genes=int(Y.shape[1]), n_p05_ols=int((p_ols < 0.05).sum()), n_p05_moderated=int((p_mod < 0.05).sum()),
           n_q05=int((q < 0.05).sum()), n_q10=int((q < 0.10).sum()), n_q20=int((q < 0.20).sum()),
           perm_null_p05_mean=float(np.mean(null)), perm_null_p05_95pct=float(np.quantile(null, 0.95)),
           empirical_FDR_p05=float(np.mean(null) / max((p_mod < 0.05).sum(), 1)), n_perm=len(null))
print(json.dumps(res, indent=1))
json.dump(res, open(os.path.join(OUT, "bulk_interaction_fdr.json"), "w"), indent=1)
pd.DataFrame({"interaction_log2": b, "p_moderated": p_mod, "q_BH": q, "p_ols": p_ols}, index=L.index).to_csv(
    os.path.join(OUT, "bulk_interaction_fdr_per_gene.tsv"), sep="\t")
