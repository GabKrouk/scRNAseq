#!/usr/bin/env python3
"""Step 26: bulk per-sample values and single-cell pseudo-bulk values (10 random
pools per state, adjusted for load/depth/cell type/other TFs) for the 35
combinatorial genes from the shared-dynamics analysis (step 24)."""
import os, sys, re
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import *
OUT = os.path.join(RES, "shared_dynamics")
sel = pd.read_csv(os.path.join(OUT, "concordant_genes_modules.tsv"), sep="\t")
sel = sel[sel.module_label.str.startswith("combination")].sort_values(["module", "score"], ascending=[True, False])
genes = sel.gene.tolist(); print(len(genes), "combinatorial concordant genes")
C = pd.read_csv("/Users/gabkrouk/Library/CloudStorage/Dropbox/Science/Papers/Emanuel New TARGET paper/Figure Will/counts.txt", sep=" ", index_col=0).drop(columns=["E2_EC.2", "E2_E.1"])
L = np.log2(C / C.sum() * 1e6 + 1).loc[genes]
bl = L.reset_index().melt(id_vars="index", var_name="sample", value_name="v").rename(columns={"index": "gene"})
bl["exp"] = bl["sample"].str[:2].map({"E1": "Bulk exp 1", "E2": "Bulk exp 2"})
bl["state"] = bl["sample"].str.extract(r"_(\w+)\.")[0].map({"N": "neither", "E": "ERF5", "C": "CRF4", "EC": "both"})
ref = bl[bl.state == "neither"].groupby(["gene", "exp"]).v.mean().rename("ref")
bl = bl.join(ref, on=["gene", "exp"]); bl["rel"] = bl.v - bl.ref
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
    B = np.linalg.lstsq(D, Y, rcond=None)[0]; adj = Y - D[:, 4:] @ B[4:]
    state = np.select([xa * xb == 1, xa == 1, xb == 1], ["both", "ERF5", "CRF4"], "neither")
    for st in ("neither", "ERF5", "CRF4", "both"):
        for k, p in enumerate(np.array_split(rng.permutation(np.where(state == st)[0]), 10)):
            for j, g in enumerate(genes):
                rows.append(dict(gene=g, exp=f"Single-cell {e}", state=st, rel=adj[p, j].mean()))
sc = pd.DataFrame(rows)
sc["rel"] = sc.rel - sc.join(sc[sc.state == "neither"].groupby(["gene", "exp"]).rel.mean().rename("r"), on=["gene", "exp"]).r
pd.concat([bl[["gene", "exp", "state", "rel"]], sc]).merge(sel[["gene", "module_name", "q", "score"]], on="gene") \
  .to_csv(os.path.join(OUT, "combinatorial_concordant_values.tsv"), sep="\t", index=False)
