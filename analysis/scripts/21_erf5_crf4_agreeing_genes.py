#!/usr/bin/env python3
"""Step 21: genes whose ERF5 x CRF4 interaction agrees between the two bulk
DoubleTARGET experiments AND the two single-cell experiments (same sign in all
four; bulk pooled p<0.05; single-cell meta p<0.1).  Exports bulk per-sample
values and single-cell pseudo-bulk values (cells of each state split at random
into 10 pools, adjusted for load/depth/cell type/other TFs) for box plots."""
import os, sys, re
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import *

OUT = os.path.join(RES, "bulk_bridge")
M = pd.read_csv(os.path.join(OUT, "interaction_sc_bulk_per_gene.tsv"), sep="\t", index_col=0)
s = np.sign
agree = (s(M.bulkE1_inter) == s(M.bulkE2_inter)) & (s(M.Exp1_inter) == s(M.Exp2_inter)) & (s(M.meta_inter) == s(M.bulk_inter))
M["tier"] = np.where(agree & (M.bulk_p_inter < .05) & (M.p_meta_inter < .05), "strict",
            np.where(agree & (M.bulk_p_inter < .05) & (M.p_meta_inter < .1), "bulk p<0.05, sc p<0.1",
            np.where(agree & (M.bulk_p_inter < .1) & (M.p_meta_inter < .1), "both p<0.1", "")))
sel = M[M.tier != ""].copy()
sel["rank"] = sel.tier.map({"strict": 0, "bulk p<0.05, sc p<0.1": 1, "both p<0.1": 2})
sel = sel.sort_values(["rank", "p_meta_inter"])
sel.to_csv(os.path.join(OUT, "agreeing_genes_ERF5xCRF4.tsv"), sep="\t")
print(sel.tier.value_counts().to_string())
genes = sel.index.tolist()

# bulk per sample (log2 CPM)
C = pd.read_csv("/Users/gabkrouk/Library/CloudStorage/Dropbox/Science/Papers/Emanuel New TARGET paper/Figure Will/counts.txt", sep=" ", index_col=0).drop(columns=["E2_EC.2", "E2_E.1"])
L = np.log2(C / C.sum() * 1e6 + 1).loc[genes]
bl = L.reset_index().melt(id_vars="index", var_name="sample", value_name="log2cpm").rename(columns={"index": "gene"})
bl["exp"] = bl["sample"].str[:2].map({"E1": "Bulk exp 1", "E2": "Bulk exp 2"})
bl["state"] = bl["sample"].str.extract(r"_(\w+)\.")[0].map({"N": "neither", "E": "ERF5", "C": "CRF4", "EC": "both"})
bl["rel"] = bl.log2cpm - bl.groupby(["gene", "exp"]).apply(lambda d: d[d.state == "neither"].log2cpm.mean()).rename("ref").reindex(pd.MultiIndex.from_frame(bl[["gene", "exp"]])).values
bl.to_csv(os.path.join(OUT, "agreeing_genes_bulk_values.tsv"), sep="\t", index=False)

# single cell pseudo-bulk
rng = np.random.default_rng(3)
ds = {e: load(e) for e in EXPS}
_, tfs = analysis_sets(ds)
CT = {e: pd.read_csv(os.path.join(RES, f"{e}_celltype_scores.tsv"), sep="\t").drop(columns=["cell_type", "cell"]) for e in EXPS}
a, b = tfs.index("ERF5"), tfs.index("CRF4")
rows = []
for e in EXPS:
    Y = response(ds[e], genes)
    X = encode(ds[e]["T"], tfs, "ge1"); xa, xb = X[:, a], X[:, b]
    S = np.column_stack([xa * (1 - xb), (1 - xa) * xb, xa * xb])
    Z = np.column_stack([covariates(ds[e], tfs, flexible_load=True).values, CT[e].values])
    D = np.column_stack([np.ones(len(xa)), S, np.delete(X, [a, b], 1), Z])
    B = np.linalg.lstsq(D, Y, rcond=None)[0]
    adj = Y - D[:, 4:] @ B[4:]
    state = np.select([xa * xb == 1, xa == 1, xb == 1], ["both", "ERF5", "CRF4"], "neither")
    for st in ("neither", "ERF5", "CRF4", "both"):
        idx = rng.permutation(np.where(state == st)[0]); pools = np.array_split(idx, 10)
        for k, p in enumerate(pools):
            for j, g in enumerate(genes):
                rows.append(dict(gene=g, exp=f"Single-cell {e}", state=st, pool=k, value=adj[p, j].mean(), n=len(p)))
sc = pd.DataFrame(rows)
ref = sc[sc.state == "neither"].groupby(["gene", "exp"]).value.mean().rename("ref")
sc = sc.join(ref, on=["gene", "exp"]); sc["rel"] = sc.value - sc.ref
sc.to_csv(os.path.join(OUT, "agreeing_genes_sc_pseudobulk.tsv"), sep="\t", index=False)
