#!/usr/bin/env python3
"""Step 52: bulk ERF5 x CRF4 DoubleTARGET - each condition vs empty vector (EV) with limma-style moderated
t-tests and BH FDR, experiment as covariate. Defines:
  ERF5-only genes   : ERF5 vs EV q<0.05, CRF4 vs EV p>0.05
  CRF4-only genes   : CRF4 vs EV q<0.05, ERF5 vs EV p>0.05
  combination-specific: ERF5+CRF4 vs EV q<0.05, and neither single TF vs EV p<0.05
  shared            : significant for a single TF and for the combination
Halving the plasmid dose weakens but does not create regulation (R3 #1), so combination-specific genes
cannot arise from the lower per-TF dose in the co-transfection."""
import os, sys, json
import numpy as np, pandas as pd
from scipy import stats, special
sys.path.insert(0, os.path.dirname(__file__))
from common import RES, bh
W = "/Users/gabkrouk/Library/CloudStorage/Dropbox/Science/Papers/Emanuel New TARGET paper/Figure Will"
OUT = os.path.join(RES, "bulk_bridge")
C = pd.read_csv(os.path.join(W, "counts.txt"), sep=" ", index_col=0).drop(columns=["E2_EC.2", "E2_E.1"])
cpm = C / C.sum() * 1e6
L = np.log2(cpm[(cpm > 1).sum(axis=1) >= 3] + 1)
s = pd.Series(L.columns); exp2 = s.str.startswith("E2").values.astype(float)
con = s.str.extract(r"_(\w+)\.")[0].values
print(pd.Series(con).value_counts().to_dict())
D = np.column_stack([np.ones(len(con)), exp2] + [(con == c).astype(float) for c in ("E", "C", "EC")])
Y = L.values.T
B = np.linalg.lstsq(D, Y, rcond=None)[0]; R = Y - D @ B
d = len(con) - D.shape[1]; s2 = (R ** 2).sum(0) / d


def inv_trigamma(x):
    y = 0.5 + 1 / x
    for _ in range(50):
        tri = special.polygamma(1, y); y = y + tri * (1 - tri / x) / special.polygamma(2, y)
    return y


e = np.log(s2) - special.digamma(d / 2) + np.log(d / 2)
v = np.var(e, ddof=1) - special.polygamma(1, d / 2)
d0 = 2 * inv_trigamma(v); s02 = np.exp(np.mean(e) + special.digamma(d0 / 2) - np.log(d0 / 2))
s2m = (d0 * s02 + d * s2) / (d0 + d)
cv = np.linalg.inv(D.T @ D)
out = pd.DataFrame(index=L.index)
for k, nm in ((2, "ERF5"), (3, "CRF4"), (4, "ERF5_CRF4")):
    t = B[k] / np.sqrt(s2m * cv[k, k]); p = 2 * stats.t.sf(np.abs(t), d + d0)
    out[f"log2FC_{nm}"] = B[k]; out[f"p_{nm}"] = p; out[f"q_{nm}"] = bh(p)
for nm, c in (("EC_vs_C", np.array([0, 0, 0, -1, 1.])), ("EC_vs_E", np.array([0, 0, -1, 0, 1.]))):
    est = c @ B; se = np.sqrt(s2m * (c @ cv @ c)); p = 2 * stats.t.sf(np.abs(est / se), d + d0)
    out[f"log2FC_{nm}"] = est; out[f"p_{nm}"] = p
qE, qC, qEC = out.q_ERF5 < .05, out.q_CRF4 < .05, out.q_ERF5_CRF4 < .05
pE, pC = out.p_ERF5 < .05, out.p_CRF4 < .05
out["class"] = np.select([qEC & ~pE & ~pC, qE & ~pC & ~qEC, qC & ~pE & ~qEC, (qE | qC) & qEC],
                         ["combination-specific", "ERF5 only", "CRF4 only", "shared single + combination"], "other")
strict = (out["class"] == "combination-specific") & (out.p_EC_vs_C < .05) & (out.p_EC_vs_E < .05)
out["combination_specific_strict"] = strict
res = {k: int(v) for k, v in out["class"].value_counts().items()}
res.update(dict(n_genes=len(out), d0=float(d0), q05_ERF5=int(qE.sum()), q05_CRF4=int(qC.sum()), q05_combination=int(qEC.sum()), combination_specific_strict=int(strict.sum()),
           comb_specific_up=int(((out["class"] == "combination-specific") & (out.log2FC_ERF5_CRF4 > 0)).sum())))
# overlap with Will's Fig 7 interaction list (1,219 genes) if available in the cluster tables
try:
    sys.path.insert(0, os.path.dirname(__file__)); from xlsx_min import read_xlsx
    for f in ("Synergy_Consensus_Clusters.xlsx", "Synergy_Non_Consensus_Clusters.xlsx"):
        rows = list(read_xlsx(os.path.join(W, f)).values())[0]; g = {r.get('A') for r in rows[1:] if r.get('A')}
        res[f"overlap_combination_specific_with_{f}"] = int(len(g & set(out.index[out["class"] == "combination-specific"])))
        res[f"n_{f}"] = len(g)
except Exception as ex:
    res["xlsx_error"] = str(ex)
print(json.dumps(res, indent=1))
json.dump(res, open(os.path.join(OUT, "bulk_condition_contrasts.json"), "w"), indent=1)
out.to_csv(os.path.join(OUT, "bulk_condition_contrasts_per_gene.tsv"), sep="\t")
