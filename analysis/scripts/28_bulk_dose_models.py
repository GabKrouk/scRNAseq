#!/usr/bin/env python3
"""Step 28: does TF dose (ERF5 / CRF4 transcript level per bulk sample) explain the
DoubleTARGET responses, in particular the ERF5 x CRF4 'interaction'?
Per gene (log2 CPM, 22 samples, E2_EC.2 and E2_E.1 removed as in Will's Rmd):
  M_fact : y ~ exp + E + C + E:C                   (factorial, dose ignored)
  M_dose : y ~ exp + dE + dC                        (continuous dose, additive)
  M_both : y ~ exp + dE + dC + E:C                  (does the interaction survive dose?)
dE / dC = log2 CPM of the ERF5 / CRF4 locus (transgene + endogenous), centred."""
import os, re, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from common import RES, bh
OUT = os.path.join(RES, "bulk_bridge")
C = pd.read_csv("/Users/gabkrouk/Library/CloudStorage/Dropbox/Science/Papers/Emanuel New TARGET paper/Figure Will/counts.txt", sep=" ", index_col=0).drop(columns=["E2_EC.2", "E2_E.1"])
cpm = C / C.sum() * 1e6
Lall = np.log2(cpm + 1)
L = Lall[(cpm > 1).sum(1) >= 3].drop(index=["AT5G47230", "AT4G27950"], errors="ignore")
s = pd.Series(L.columns)
con = s.str.extract(r"_(\w+)\.")[0].values
exp2 = s.str.startswith("E2").astype(float).values
E = np.isin(con, ["E", "EC"]).astype(float); Cc = np.isin(con, ["C", "EC"]).astype(float)
dE = Lall.loc["AT5G47230"].values; dC = Lall.loc["AT4G27950"].values
dose = pd.DataFrame({"sample": L.columns, "construct": con, "ERF5_log2CPM": dE, "CRF4_log2CPM": dC})
dose.to_csv(os.path.join(OUT, "bulk_TF_dose_per_sample.tsv"), sep="\t", index=False)
print(dose.groupby("construct")[["ERF5_log2CPM", "CRF4_log2CPM"]].mean().round(2))
dEc, dCc = dE - dE.mean(), dC - dC.mean()
Y = L.values.T; n = Y.shape[0]


def fit(D, k=None):
    B, *_ = np.linalg.lstsq(D, Y, rcond=None); R = Y - D @ B
    rss = (R ** 2).sum(0); df = n - D.shape[1]
    out = dict(rss=rss, df=df, B=B, aic=n * np.log(rss / n) + 2 * D.shape[1])
    if k is not None:
        se = np.sqrt(rss / df * np.linalg.inv(D.T @ D)[k, k]); t = B[k] / se
        out["t"] = t; out["p"] = 2 * stats.t.sf(np.abs(t), df)
    return out


one = np.ones(n)
fact = fit(np.column_stack([one, exp2, E, Cc, E * Cc]), 4)
dosem = fit(np.column_stack([one, exp2, dEc, dCc]))
both = fit(np.column_stack([one, exp2, dEc, dCc, E * Cc]), 4)
both_noint = dosem
# does dose explain the factorial pattern? compare M_fact vs M_dose by AIC (same #params minus 1)
res = pd.DataFrame(index=L.index)
res["inter_fact"] = fact["B"][4]; res["p_inter_fact"] = fact["p"]
res["inter_after_dose"] = both["B"][4]; res["p_inter_after_dose"] = both["p"]
res["aic_fact"] = fact["aic"]; res["aic_dose"] = dosem["aic"]
res["dose_better"] = res.aic_dose < res.aic_fact - 2
res["fact_better"] = res.aic_fact < res.aic_dose - 2
res.to_csv(os.path.join(OUT, "bulk_dose_models_per_gene.tsv"), sep="\t")
print("genes:", len(res))
print("interaction p<0.05 factorial: %d (%.1f%%) | after dose covariates: %d (%.1f%%) | expected by chance 5%%" % (
    (res.p_inter_fact < .05).sum(), 100 * (res.p_inter_fact < .05).mean(), (res.p_inter_after_dose < .05).sum(), 100 * (res.p_inter_after_dose < .05).mean()))
print("corr interaction estimate before/after dose: %.3f" % np.corrcoef(res.inter_fact, res.inter_after_dose)[0, 1])
print("AIC: dose model better (dAIC>2): %d | factorial better: %d | undecided: %d" % (res.dose_better.sum(), res.fact_better.sum(), (~res.dose_better & ~res.fact_better).sum()))
# among genes that actually respond (any effect), restricted view
resp = fit(np.column_stack([one, exp2]))
F = ((resp["rss"] - fact["rss"]) / 3) / (fact["rss"] / fact["df"]); pF = stats.f.sf(F, 3, fact["df"])
res["p_any_response"] = pF; res["q_any_response"] = bh(pF)
r = res[res.q_any_response < 0.05]
print("responsive genes (any construct effect, q<0.05):", len(r))
print("  among them: dose model better %d | factorial better %d | interaction p<0.05 factorial %d -> after dose %d" % (
    r.dose_better.sum(), r.fact_better.sum(), (r.p_inter_fact < .05).sum(), (r.p_inter_after_dose < .05).sum()))
# shared-dynamics genes
sd = pd.read_csv(os.path.join(RES, "shared_dynamics", "concordant_genes_modules.tsv"), sep="\t")
for lab, g in (("333 concordant", sd), ("35 combinatorial", sd[sd.module_label.str.startswith("combination")])):
    x = res.reindex(g.gene).dropna()
    print(f"{lab}: n={len(x)} | dose model better {x.dose_better.sum()} | factorial better {x.fact_better.sum()} | interaction p<0.05: factorial {(x.p_inter_fact<.05).sum()} -> after dose {(x.p_inter_after_dose<.05).sum()}")
res.to_csv(os.path.join(OUT, "bulk_dose_models_per_gene.tsv"), sep="\t")
# CRF4 dose-response among strong CRF4 targets: slope of target response on CRF4 level within CRF4-containing samples
crf = res[(res.q_any_response < .05)].index
