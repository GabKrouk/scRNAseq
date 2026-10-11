#!/usr/bin/env python3
"""Step 48: tidy tables for the revised paper figures (raw counts, QC cells).
TF landscape: detection per TF, number of TFs per cell, pairwise co-detection
(observed vs expected under independence), per experiment."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import *
OUT = os.path.join(RES, "paper"); os.makedirs(OUT, exist_ok=True)
ds = {e: load(e) for e in EXPS}
genes, tfs = analysis_sets(ds)
det, ntf, co = [], [], []
for e in EXPS:
    B = (ds[e]["T"][tfs].values >= 1)
    n = B.sum(1)
    for k in range(len(tfs) + 1):
        ntf.append(dict(exp=e, n_tfs=k, n_cells=int((n == k).sum()), frac=float((n == k).mean())))
    for j, t in enumerate(tfs):
        det.append(dict(exp=e, tf=t, n_pos=int(B[:, j].sum()), frac=float(B[:, j].mean())))
    p = B.mean(0)
    for a in range(len(tfs)):
        for b in range(len(tfs)):
            o = float((B[:, a] & B[:, b]).mean())
            co.append(dict(exp=e, a=tfs[a], b=tfs[b], obs=o, exp_indep=float(p[a] * p[b]) if a != b else float(p[a]),
                           log2_oe=float(np.log2(o / (p[a] * p[b]))) if a != b and o > 0 else np.nan))
    print(e, "cells", len(n), "median TFs/cell", np.median(n), "no TF", (n == 0).mean().round(3), "mean", n.mean().round(2))
pd.DataFrame(det).to_csv(os.path.join(OUT, "tf_detection.tsv"), sep="\t", index=False)
pd.DataFrame(ntf).to_csv(os.path.join(OUT, "tfs_per_cell.tsv"), sep="\t", index=False)
pd.DataFrame(co).to_csv(os.path.join(OUT, "tf_codetection.tsv"), sep="\t", index=False)
