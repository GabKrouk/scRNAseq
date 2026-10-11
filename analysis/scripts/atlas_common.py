"""Loader for the Shahan et al. (2022) WT root atlas (raw UMIs from step 65) and the fixed-effect machinery
of steps 66-69. Default stratum = sample x cell type x developmental time group; covariates = log depth,
(log depth)^2, log genes detected, % mitochondrial reads. make_tilde() builds the same machinery for any
cell subset and any stratum definition."""
import os
import numpy as np, pandas as pd
from scipy import sparse
from common import TF_AGI

ATLAS = os.path.expanduser("~/scdata/Shahan2022")
AGI = dict(TF_AGI); AGI["HHO3"] = "AT1G25550"      # endogenous HHO3 (NIGT1.1) locus
STRATA_FULL = ["orig.ident", "celltype.anno", "consensus.time.group"]


def make_tilde(meta, tot, ngen, mask, keys):
    """Residualiser on stratum fixed effects (keys) + covariates, for the cells in mask.
    Returns (tilde, DF); tilde takes arrays over the masked cells (n_mask x k)."""
    m = meta[mask]
    strat = m[keys[0]].astype(str).values if keys else np.zeros(mask.sum(), dtype=str)
    for k in keys[1:]:
        strat = strat + "|" + m[k].astype(str).values
    codes, levels = pd.factorize(strat); K = len(levels); n = mask.sum()
    cnt = np.bincount(codes, minlength=K).astype(np.float64)
    S = sparse.csr_matrix((np.ones(n), (np.arange(n), codes)), shape=(n, K))

    def demean(X):
        X = np.asarray(X, np.float64)
        return X - ((S.T @ X) / cnt[:, None])[codes]

    lt = np.log(tot[mask])
    C = demean(np.column_stack([lt, lt ** 2, np.log(ngen[mask]), m["percent.mt"].astype(float).values]))
    Cp = np.linalg.pinv(C)

    def tilde(X):
        X = demean(X)
        return X - C @ (Cp @ X)

    return tilde, n - K - C.shape[1]


def load_atlas():
    D = os.path.join(ATLAS, "derived")
    ng, nc = [int(v) for v in open(os.path.join(D, "dim.txt")).read().split()]
    M = sparse.csc_matrix((np.fromfile(os.path.join(D, "counts_x.int32"), np.int32).astype(np.float32),
                           np.fromfile(os.path.join(D, "counts_i.int32"), np.int32),
                           np.fromfile(os.path.join(D, "counts_p.int32"), np.int32)), shape=(ng, nc))
    genes = pd.Index(open(os.path.join(D, "genes.txt")).read().split())
    meta = pd.read_csv(os.path.join(ATLAS, "metadata", "Shahan2022_atlas_cell_metadata.tsv.gz"), sep="\t",
                       low_memory=False).set_index("cell_barcode_id")
    cc = pd.read_csv(os.path.join(D, "cells_counts.tsv"), sep="\t")
    same = ((cc.cell.str[:16].values == meta.index.str[:16].values)
            & (cc.nCount_spliced.values == meta.nCount_spliced_RNA.astype(int).values)
            & (cc.nCount_unspliced.values == meta.nCount_unspliced_RNA.astype(int).values))
    assert same.all()
    tot = np.asarray(M.sum(0)).ravel().astype(np.float64)
    ngen = np.diff(M.indptr).astype(np.float64)
    tilde, DF = make_tilde(meta, tot, ngen, np.ones(nc, bool), STRATA_FULL)
    return dict(M=M, genes=genes, meta=meta, tot=tot, ngen=ngen, n=nc, tilde=tilde, DF=DF)
