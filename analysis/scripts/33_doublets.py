#!/usr/bin/env python3
"""Step 33: doublet scores (Scrublet-like, implemented in numpy): simulate doublets by
summing the raw counts of random cell pairs of the same library, embed real and
simulated profiles in a PCA fitted on real cells, and score each real cell by the
fraction of simulated doublets among its k nearest neighbours. Cells are called
doublets at the expected 10x multiplet rate (0.8% per 1,000 recovered cells, per library)."""
import os, sys
import numpy as np, pandas as pd, scipy.io as sio
from scipy import sparse
sys.path.insert(0, os.path.dirname(__file__))
from common import RES
TOP = "/Users/gabkrouk/Library/CloudStorage/Dropbox/Science/projects/analyse_gab_scTARGET"
OUT = os.path.join(RES, "qc_raw"); rng = np.random.default_rng(17)
for e in ("Exp1", "Exp2"):
    qc = pd.read_csv(os.path.join(OUT, f"{e}_cell_qc.tsv.gz"), sep="\t")
    M = sio.mmread(os.path.join(TOP, "raw_data", f"{e}_matrix.mtx")).tocsc()
    feat = pd.read_csv(os.path.join(TOP, f"{e}_features.tsv.gz"), sep="\t", header=None, names=["id", "name", "type"])
    keepf = feat.id.str.match(r"^AT[1-5]G\d{5}").values
    X = M[keepf].T.tocsr().astype(np.float32)                 # cells x genes
    scores = np.zeros(X.shape[0]); calls = np.zeros(X.shape[0], bool)
    for lib, idx in qc.groupby("library").groups.items():
        idx = np.asarray(idx); Xl = X[idx]; n = len(idx)
        nsim = 2 * n
        a, b = rng.integers(0, n, nsim), rng.integers(0, n, nsim)
        Xs = Xl[a] + Xl[b]
        def norm(Z):
            tot = np.asarray(Z.sum(1)).ravel()
            return sparse.diags(1e4 / tot) @ Z
        Ln = norm(Xl); Ls = norm(Xs)
        Ln.data = np.log1p(Ln.data); Ls.data = np.log1p(Ls.data)
        mu = np.asarray(Ln.mean(0)).ravel(); var = np.asarray(Ln.multiply(Ln).mean(0)).ravel() - mu ** 2
        hv = np.argsort(var)[-2000:]
        A = Ln[:, hv].toarray(); S = Ls[:, hv].toarray()
        m, sd = A.mean(0), A.std(0) + 1e-6
        A = (A - m) / sd; S = (S - m) / sd
        U, Sv, Vt = np.linalg.svd(A, full_matrices=False)
        P = A @ Vt[:30].T; Q = S @ Vt[:30].T
        allp = np.vstack([P, Q]); lab = np.r_[np.zeros(n), np.ones(nsim)]
        k = max(10, int(round(0.5 * np.sqrt(n))))
        d2 = (P ** 2).sum(1)[:, None] + (allp ** 2).sum(1)[None, :] - 2 * P @ allp.T
        d2[np.arange(n), np.arange(n)] = np.inf
        nn = np.argpartition(d2, k, axis=1)[:, :k]
        frac = lab[nn].mean(1)
        rate = 0.008 * n / 1000
        ncall = int(np.ceil(rate * n))
        order = np.lexsort((rng.random(n), -frac))[:ncall]          # top scores, ties broken at random
        call = np.zeros(n, bool); call[order] = True; call &= frac > 0
        scores[idx] = frac; calls[idx] = call
        print(e, "lib", lib, "cells", n, "expected rate %.1f%%" % (100 * rate), "called", int(call.sum()),
              "| median score %.2f, score of called >= %.2f" % (np.median(frac), frac[call].min() if call.any() else np.nan))
    qc["doublet_score"] = scores; qc["doublet"] = calls
    qc.to_csv(os.path.join(OUT, f"{e}_cell_qc.tsv.gz"), sep="\t", index=False)
