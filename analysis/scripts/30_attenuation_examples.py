#!/usr/bin/env python3
"""Step 30: example genes for the dose argument (Reviewer 3 #1).
CRF4-induced genes whose induction is attenuated by ERF5: in bulk this is
explained by the halved CRF4 dose (dose model wins), but in single cells, where
CRF4 dose is the same in CRF4-only and ERF5+CRF4 cells, the attenuation persists.
Ranks module 'CRF4-driven induction' genes by the single-cell attenuation
(both - CRF4 only, HC3 t) in both experiments; exports values for box plots."""
import os, sys, re
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from common import *
OUT = os.path.join(RES, "bulk_bridge")
sd = pd.read_csv(os.path.join(RES, "shared_dynamics", "concordant_genes_modules.tsv"), sep="\t")
dm = pd.read_csv(os.path.join(OUT, "bulk_dose_models_per_gene.tsv"), sep="\t", index_col=0)
cand = sd[sd.module_label == "CRF4-driven induction"].gene.tolist()
ds = {e: load(e) for e in EXPS}
_, tfs = analysis_sets(ds)
CT = {e: pd.read_csv(os.path.join(RES, f"{e}_celltype_scores.tsv"), sep="\t").drop(columns=["cell_type", "cell"]) for e in EXPS}
a, b = tfs.index("ERF5"), tfs.index("CRF4")
rng = np.random.default_rng(4)
stat, cells, dose_sc = {}, {}, []
for e in EXPS:
    Y = response(ds[e], cand)
    X = encode(ds[e]["T"], tfs, "ge1"); xa, xb = X[:, a], X[:, b]
    S = np.column_stack([xa * (1 - xb), (1 - xa) * xb, xa * xb])          # ERF5 only, CRF4 only, both
    Z = np.column_stack([covariates(ds[e], tfs, flexible_load=True).values, CT[e].values])
    D = np.column_stack([np.ones(len(xa)), S, np.delete(X, [a, b], 1), Z])
    A = np.linalg.pinv(D.T @ D) @ D.T; B = A @ Y; E = Y - D @ B
    h = np.einsum("ij,ji->i", D, A); W = (E / (1 - h)[:, None]) ** 2
    c = A[3] - A[2]                                                       # both - CRF4 only
    est = B[3] - B[2]; se = np.sqrt((c ** 2) @ W)
    stat[e] = pd.DataFrame({"atten": est, "t": est / se, "crf4_effect": B[2]}, index=cand)
    adj = Y - D[:, 4:] @ B[4:]
    state = np.select([xa * xb == 1, xa == 1, xb == 1], ["both", "ERF5", "CRF4"], "neither")
    cells[e] = (adj, state)
    crf = ds[e]["T"].CRF4.values
    for st in ("CRF4", "both"):
        m = state == st
        dose_sc.append(pd.DataFrame({"dataset": f"Single-cell {e}", "state": st, "dose": np.log1p(crf[m])}))
res = pd.DataFrame(index=cand)
for e in EXPS:
    res[f"atten_{e}"] = stat[e].atten; res[f"t_{e}"] = stat[e].t
res["bulk_dose_better"] = dm.reindex(cand).dose_better.values
res["min_neg_t"] = -res[[f"t_{e}" for e in EXPS]].max(1)                 # attenuation = negative t in both
sel = res[(res[[f"t_{e}" for e in EXPS]] < -1.96).all(1) & res.bulk_dose_better].sort_values("min_neg_t", ascending=False)
print("CRF4-induced genes:", len(cand), "| attenuated (t<-1.96) in both sc exps:", ((res[[f't_{e}' for e in EXPS]] < -1.96).all(1)).sum(),
      "| of which bulk dose model wins:", len(sel))
print(sel.head(10).round(3).to_string())
genes = res[(res[[f"t_{e}" for e in EXPS]] < -1.96).all(1)].sort_values("min_neg_t", ascending=False).index.tolist()
res.to_csv(os.path.join(OUT, "attenuation_all_CRF4_genes.tsv"), sep="\t"); sel.to_csv(os.path.join(OUT, "attenuation_examples_ranked.tsv"), sep="\t")
# values for plots
rows = []
for e in EXPS:
    adj, state = cells[e]
    for st in ("neither", "ERF5", "CRF4", "both"):
        idx = rng.permutation(np.where(state == st)[0])
        for k, p in enumerate(np.array_split(idx, 10)):
            for g in genes:
                rows.append(dict(gene=g, exp=f"Single-cell {e}", state=st, v=adj[p, cand.index(g)].mean()))
sc = pd.DataFrame(rows)
sc["rel"] = sc.v - sc.join(sc[sc.state == "neither"].groupby(["gene", "exp"]).v.mean().rename("r"), on=["gene", "exp"]).r
Cn = pd.read_csv("/Users/gabkrouk/Library/CloudStorage/Dropbox/Science/Papers/Emanuel New TARGET paper/Figure Will/counts.txt", sep=" ", index_col=0).drop(columns=["E2_EC.2", "E2_E.1"])
Lb = np.log2(Cn / Cn.sum() * 1e6 + 1)
bl = Lb.loc[genes].reset_index().melt(id_vars="index", var_name="sample", value_name="v").rename(columns={"index": "gene"})
bl["exp"] = bl["sample"].str[:2].map({"E1": "Bulk exp 1", "E2": "Bulk exp 2"})
bl["state"] = bl["sample"].str.extract(r"_(\w+)\.")[0].map({"N": "neither", "E": "ERF5", "C": "CRF4", "EC": "both"})
bl["rel"] = bl.v - bl.join(bl[bl.state == "neither"].groupby(["gene", "exp"]).v.mean().rename("r"), on=["gene", "exp"]).r
pd.concat([bl[["gene", "exp", "state", "rel"]], sc[["gene", "exp", "state", "rel"]]]).to_csv(os.path.join(OUT, "attenuation_examples_values.tsv"), sep="\t", index=False)
db = Lb.loc["AT4G27950"].reset_index(); db.columns = ["sample", "dose"]
db["dataset"] = db["sample"].str[:2].map({"E1": "Bulk exp 1", "E2": "Bulk exp 2"})
db["state"] = db["sample"].str.extract(r"_(\w+)\.")[0].map({"N": "neither", "E": "ERF5", "C": "CRF4", "EC": "both"})
db = db[db.state.isin(["CRF4", "both"])]
pd.concat([db[["dataset", "state", "dose"]], pd.concat(dose_sc)]).to_csv(os.path.join(OUT, "attenuation_examples_CRF4dose.tsv"), sep="\t", index=False)
for e in EXPS:
    d = pd.concat(dose_sc); d = d[d.dataset == f"Single-cell {e}"]
    print(e, "CRF4 dose (log1p reads) CRF4-only vs both: %.2f vs %.2f, MW p=%.2g" % (d[d.state == "CRF4"].dose.mean(), d[d.state == "both"].dose.mean(),
          stats.mannwhitneyu(d[d.state == "CRF4"].dose, d[d.state == "both"].dose).pvalue))
