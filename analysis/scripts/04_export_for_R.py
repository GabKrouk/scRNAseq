#!/usr/bin/env python3
"""Step 4: export matrices for the CART analyses in R (rpart).
For each experiment (cells x genes, column-major float32 / int32 raw files):
  *_Yres.f32   log1p-CP10k expression residualised on covariates Z (depth, genes
               detected, RFP, total GR-TF reads incl. non-linear terms, library,
               marker-based cell-type scores)  -> TF-attributable signal
  *_Yraw.i32   raw integer counts (for the CART-CELLS-faithful classification mode)
  *_Ylog.f32   log1p-CP10k (no covariate removal)
  *_TF.i32     GR-TF read counts for the analysed TFs
  *_Z.f32      covariates (for partialling predictions)
"""
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import *

OUTR = os.path.join(DATA, "forR"); os.makedirs(OUTR, exist_ok=True)
ds = {e: load(e) for e in EXPS}
genes, tfs = analysis_sets(ds)
meta = dict(genes=genes, tfs=tfs, n={})
for e in EXPS:
    d = ds[e]
    Y = response(d, genes)
    ct = pd.read_csv(os.path.join(RES, f"{e}_celltype_scores.tsv"), sep="\t").drop(columns=["cell_type", "cell"])
    Z = pd.concat([covariates(d, tfs, flexible_load=True), ct], axis=1)
    D = np.column_stack([np.ones(len(Y)), Z.values])
    Yres = Y - D @ np.linalg.lstsq(D, Y, rcond=None)[0]
    idx = pd.Index(d["genes"]).get_indexer(genes)
    Yraw = d["G"][idx, :].toarray().T.astype(np.int32)
    TF = d["T"][tfs].values.astype(np.int32)
    np.asfortranarray(Yres.astype(np.float32)).T.tofile(os.path.join(OUTR, f"{e}_Yres.f32"))
    np.asfortranarray(Y.astype(np.float32)).T.tofile(os.path.join(OUTR, f"{e}_Ylog.f32"))
    np.asfortranarray(Yraw).T.tofile(os.path.join(OUTR, f"{e}_Yraw.i32"))
    np.asfortranarray(TF).T.tofile(os.path.join(OUTR, f"{e}_TF.i32"))
    np.asfortranarray(Z.values.astype(np.float32)).T.tofile(os.path.join(OUTR, f"{e}_Z.f32"))
    meta["n"][e] = dict(cells=int(Y.shape[0]), z_cols=Z.columns.tolist())
    print(e, Y.shape, TF.shape, Z.shape)
pd.Series(genes).to_csv(os.path.join(OUTR, "genes.txt"), index=False, header=False)
pd.Series(tfs).to_csv(os.path.join(OUTR, "tfs.txt"), index=False, header=False)
json.dump(meta["n"], open(os.path.join(OUTR, "dims.json"), "w"), indent=1)
