"""Shared loaders and helpers for the scTARGET-seq reproducibility analysis."""
import os
import numpy as np, pandas as pd
from scipy import sparse

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data")
RES = os.path.join(BASE, "results")
FIG = os.path.join(BASE, "figures")
EXPS = ["Exp1", "Exp2"]

# AGI of the endogenous locus of each transfected TF (excluded as targets so
# transgene-derived reads cannot create self "targets"). Checked against the
# data in 02_explore.py (TF-row vs locus correlation).
TF_AGI = {
    "BEE2": "AT4G36540", "CDF1": "AT5G62430", "bZIP3": "AT5G15830",
    "COL5": "AT5G57660", "CRF4": "AT4G27950", "ERF5": "AT5G47230",
    "FBH4": "AT2G42280", "GATA17": "AT3G16870", "GATA17L": "AT4G16141",
    "HAT22": "AT4G37790", "HB6": "AT2G22430", "HHO2": "AT1G68670",
    "HHO3": "AT2G03500", "HSFB2A": "AT5G62020", "HYH": "AT3G17609",
    "LBD37": "AT5G67420", "LBD38": "AT3G49940", "MYB61": "AT1G09540",
    "NAC4": "AT5G07680", "NAP": "AT1G69490", "RAV1": "AT1G13260",
    "TCP23": "AT1G35560", "TGA1": "AT5G65210", "TGA4": "AT5G10030",
    "VRN1": "AT3G18990", "WRKY18": "AT4G31800", "WRKY54": "AT2G40750",
}


def load(exp):
    G = sparse.load_npz(os.path.join(DATA, f"{exp}_genes.npz")).tocsc()
    genes = pd.read_csv(os.path.join(DATA, f"{exp}_genes.txt"), header=None)[0].tolist()
    T = np.load(os.path.join(DATA, f"{exp}_tf.npy"))
    tfs = pd.read_csv(os.path.join(DATA, f"{exp}_tfnames.txt"), header=None)[0].tolist()
    rfp = np.load(os.path.join(DATA, f"{exp}_rfp.npy"))
    cells = pd.read_csv(os.path.join(DATA, f"{exp}_cells.tsv"), sep="\t", dtype=str)
    return dict(G=G, genes=genes, T=pd.DataFrame(T.T, columns=tfs), rfp=rfp, cells=cells)


def bh(p):
    """Benjamini-Hochberg q-values for an array of any shape (NaN-safe)."""
    p = np.asarray(p, float)
    flat = p.ravel()
    ok = ~np.isnan(flat)
    q = np.full(flat.shape, np.nan)
    pv = flat[ok]
    n = pv.size
    o = np.argsort(pv)
    r = pv[o] * n / np.arange(1, n + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    tmp = np.empty(n); tmp[o] = np.minimum(r, 1)
    q[ok] = tmp
    return q.reshape(p.shape)


# ---------------------------------------------------------------- analysis set
MIN_DETECT = 0.05          # gene detected in >=5% of cells in BOTH experiments
MIN_TF_CELLS = 100         # TF must be detected (>=1 read) in >=100 cells in BOTH


def analysis_sets(ds):
    """Genes and TFs kept for all analyses (identical across experiments)."""
    fr = {e: pd.Series((d["G"] > 0).mean(1).A1, index=d["genes"]) for e, d in ds.items()}
    common = fr["Exp1"].index.intersection(fr["Exp2"].index)
    keep = common[(fr["Exp1"][common] >= MIN_DETECT) & (fr["Exp2"][common] >= MIN_DETECT)]
    tf_loci = set(TF_AGI.values())
    keep = [g for g in keep if g not in tf_loci]
    tfs = [t for t in ds["Exp1"]["T"].columns if t in ds["Exp2"]["T"].columns
           and all((ds[e]["T"][t] >= 1).sum() >= MIN_TF_CELLS for e in EXPS)]
    return keep, tfs


def response(d, genes):
    """log1p(CP10k) over endogenous genes (TF rows and RFP excluded from depth)."""
    G = d["G"]
    idx = pd.Index(d["genes"]).get_indexer(genes)
    tot = G.sum(0).A1
    Y = G[idx, :].toarray().T.astype(np.float64)
    Y = np.log1p(Y / tot[:, None] * 1e4)
    return Y


def covariates(d, tfs, include_library=True, extra=None, flexible_load=False):
    """Technical / load covariates: depth, genes detected, RFP, total GR-TF reads."""
    G, T = d["G"], d["T"]
    Z = pd.DataFrame({
        "log_umi": np.log(G.sum(0).A1),
        "log_ngenes": np.log((G > 0).sum(0).A1),
        "log_rfp": np.log1p(d["rfp"]),
        "log_tf_total": np.log1p(T.sum(1).values),
    })
    if flexible_load:
        # non-linear load terms: without them the TF x TF product terms pick up
        # non-linear load effects (interaction test inflated under permutation,
        # step 5). Polynomials of n_TF are NOT used: n_TF is a function of the
        # TF design itself and absorbs true TF signal without improving the null.
        for c in ("log_umi", "log_rfp", "log_tf_total"):
            v = (Z[c] - Z[c].mean()) / Z[c].std()
            Z[c + "^2"] = v ** 2
            Z[c + "^3"] = v ** 3
        Z["logTF_x_logRFP"] = ((Z["log_tf_total"] - Z["log_tf_total"].mean()) *
                               (Z["log_rfp"] - Z["log_rfp"].mean()))
    if include_library and d["cells"]["library"].nunique() > 1:
        L = pd.get_dummies(d["cells"]["library"], prefix="lib", drop_first=True).astype(float)
        Z = pd.concat([Z, L.reset_index(drop=True)], axis=1)
    if extra is not None:
        Z = pd.concat([Z, extra.reset_index(drop=True)], axis=1)
    return Z


def encode(T, tfs, scheme):
    """TF design matrix. scheme: 'ge1','ge2','ge3','ge5' (binary), 'log' (log1p
    counts), 'share' (log share of the cell's GR-TF reads, pseudo-count 0.5)."""
    M = T[tfs].values.astype(float)
    if scheme.startswith("ge"):
        return (M >= int(scheme[2:])).astype(float)
    if scheme == "log":
        return np.log1p(M)
    if scheme == "share":
        tot = T.values.sum(1, keepdims=True)
        return np.log((M + 0.5) / (tot + 0.5 * T.shape[1]))
    raise ValueError(scheme)


def ols(X, Z, Y, robust=True):
    """Vectorised OLS of every gene (columns of Y) on [1, X, Z].
    Returns beta, se, t for the X columns, the residuals and the design."""
    from scipy import stats
    n = Y.shape[0]
    D = np.column_stack([np.ones(n), X, Z]) if Z is not None and Z.shape[1] else np.column_stack([np.ones(n), X])
    XtX_inv = np.linalg.pinv(D.T @ D)
    A = XtX_inv @ D.T                       # p x n
    B = A @ Y                               # p x G
    E = Y - D @ B
    p = D.shape[1]
    k = X.shape[1]
    if robust:                              # HC3 sandwich (leverage-corrected)
        h = np.einsum("ij,ji->i", D, A)     # leverages diag(D (D'D)^-1 D')
        V = (A[1:k + 1] ** 2) @ ((E / (1 - h)[:, None]) ** 2)
    else:
        s2 = (E ** 2).sum(0) / (n - p)
        V = np.outer(np.diag(XtX_inv)[1:k + 1], s2)
    se = np.sqrt(V)
    t = B[1:k + 1] / se
    pval = 2 * stats.t.sf(np.abs(t), n - p)
    return dict(beta=B[1:k + 1], se=se, t=t, p=pval, resid=E, design=D, df=n - p,
                lev=h if robust else None)


def strata_permutation(d, rng, nbins=10):
    """Permutation of cells' TF profiles within strata of library x total-GR-TF
    -read decile (keeps the load structure, breaks TF identity -> expression)."""
    tot = d["T"].sum(1).values
    q = pd.qcut(np.log1p(tot) + rng.uniform(0, 1e-6, tot.size), nbins, labels=False, duplicates="drop")
    key = d["cells"]["library"].astype(str).values + "_" + q.astype(str)
    perm = np.arange(tot.size)
    for k in np.unique(key):
        idx = np.where(key == k)[0]
        perm[idx] = rng.permutation(idx)
    return perm
