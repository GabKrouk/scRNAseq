#!/usr/bin/env python3
"""Step 1 (raw-count version): build the analysis data set from the raw CellRanger
7.1.0 filtered matrices (dataverse_files.zip: Exp1 2,316 cells, Exp2 6,572 cells,
reference = Col-0 genome + GR-TF fusion junctions + RFP + GFP).

Cell QC (per 10x library):
  - low UMI or low genes: < median - 3 MAD on log10 scale
  - mitochondrial (ATMG + atp4) > 5% or plastid (ATCG) > 5% of UMIs
  - doublets: Scrublet-like simulation, called at the expected 10x multiplet rate
    (0.8% per 1,000 recovered cells)
Only QC-passing cells are written to ../data/ in the format used by all later steps
(genes x cells sparse counts, GR-TF junction counts, RFP counts, cell table).
"""
import os, json
import numpy as np, pandas as pd, scipy.io as sio
from scipy import sparse

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = "/Users/gabkrouk/Library/CloudStorage/Dropbox/Science/projects/analyse_gab_scTARGET/dataverse"
OUT = os.path.join(BASE, "data"); QC = os.path.join(BASE, "results", "qc"); os.makedirs(QC, exist_ok=True)
rng = np.random.default_rng(17)


def mad_low(x):
    lx = np.log10(x); med = np.median(lx); mad = 1.4826 * np.median(np.abs(lx - med))
    return lx < med - 3 * mad


def doublet_scores(Xl):
    """Xl: cells x genes raw counts (csr). Returns doublet score per cell."""
    n = Xl.shape[0]; nsim = 2 * n
    a, b = rng.integers(0, n, nsim), rng.integers(0, n, nsim)
    Xs = Xl[a] + Xl[b]
    def lnorm(Z):
        Z = sparse.diags(1e4 / np.asarray(Z.sum(1)).ravel()) @ Z
        Z = Z.tocsr(); Z.data = np.log1p(Z.data); return Z
    Ln, Ls = lnorm(Xl), lnorm(Xs)
    mu = np.asarray(Ln.mean(0)).ravel(); var = np.asarray(Ln.multiply(Ln).mean(0)).ravel() - mu ** 2
    hv = np.argsort(var)[-2000:]
    A, S = Ln[:, hv].toarray(), Ls[:, hv].toarray()
    m, sd = A.mean(0), A.std(0) + 1e-6
    A, S = (A - m) / sd, (S - m) / sd
    Vt = np.linalg.svd(A, full_matrices=False)[2][:30]
    P, Q = A @ Vt.T, S @ Vt.T
    allp = np.vstack([P, Q]); lab = np.r_[np.zeros(n), np.ones(nsim)]
    k = max(10, int(round(0.5 * np.sqrt(n))))
    d2 = (P ** 2).sum(1)[:, None] + (allp ** 2).sum(1)[None, :] - 2 * P @ allp.T
    d2[np.arange(n), np.arange(n)] = np.inf
    return lab[np.argpartition(d2, k, axis=1)[:, :k]].mean(1)


summary = {}
for e in ("Exp1", "Exp2"):
    M = sio.mmread(os.path.join(SRC, f"{e}_matrix.mtx.gz")).tocsr()
    feat = pd.read_csv(os.path.join(SRC, f"{e}_features.tsv.gz"), sep="\t", header=None, names=["id", "name", "type"])
    bc = pd.read_csv(os.path.join(SRC, f"{e}_barcodes.tsv.gz"), header=None)[0]
    isgene = (feat.id.str.match(r"^AT[1-5CM]G\d{5}") | (feat.name == "atp4")).values
    mito = (feat.id.str.startswith("ATMG") | (feat.name == "atp4")).values
    plast = feat.id.str.startswith("ATCG").values
    # look up transgene features by feature ID (not by name: gene symbols such as ERF5 also exist)
    tf_idx = [i for i in np.where(~isgene)[0] if feat.id[i] not in ("RFP", "GFP")]
    tf_names = feat.id[tf_idx].tolist()
    G = M[isgene]
    umi = np.asarray(G.sum(0)).ravel(); ngen = np.asarray((G > 0).sum(0)).ravel()
    q = pd.DataFrame({"cell": bc.values, "library": bc.str.split("-").str[-1].values, "umi": umi, "genes": ngen,
                      "pct_mito": 100 * np.asarray(M[mito].sum(0)).ravel() / umi,
                      "pct_plastid": 100 * np.asarray(M[plast].sum(0)).ravel() / umi,
                      "rfp": np.asarray(M[np.where(feat.id.values == "RFP")[0][0]].todense()).ravel(),
                      "gfp": np.asarray(M[np.where(feat.id.values == "GFP")[0][0]].todense()).ravel()})
    T = np.asarray(M[tf_idx].todense())
    q["tf_reads"] = T.sum(0); q["n_tf"] = (T > 0).sum(0)
    q["low_umi"] = False; q["low_genes"] = False; q["doublet_score"] = np.nan; q["doublet"] = False
    Xc = M[isgene & ~mito & ~plast].T.tocsr().astype(np.float32)
    for lib, idx in q.groupby("library").groups.items():
        idx = np.asarray(idx)
        q.loc[idx, "low_umi"] = mad_low(q.umi.values[idx]); q.loc[idx, "low_genes"] = mad_low(q.genes.values[idx])
        s = doublet_scores(Xc[idx]); n = len(idx); ncall = int(np.ceil(0.008 * n / 1000 * n))
        call = np.zeros(n, bool); call[np.lexsort((rng.random(n), -s))[:ncall]] = True; call &= s > 0
        q.loc[idx, "doublet_score"] = s; q.loc[idx, "doublet"] = call
    q["high_mito"] = q.pct_mito > 5; q["high_plastid"] = q.pct_plastid > 5
    q["pass_qc"] = ~(q.low_umi | q.low_genes | q.high_mito | q.high_plastid | q.doublet)
    # cells used in the paper's normalised matrix (for reporting only)
    paper = pd.read_csv(os.path.join(BASE, "..", "reproducibility", "data", f"{e}_cells.tsv"), sep="\t").cell.str.replace(".", "-", regex=False)
    q["in_paper_matrix"] = q.cell.isin(set(paper))
    q.to_csv(os.path.join(QC, f"{e}_cell_qc.tsv.gz"), sep="\t", index=False)
    k = q.pass_qc.values
    sparse.save_npz(os.path.join(OUT, f"{e}_genes.npz"), sparse.csr_matrix(G[:, k].astype(np.float32)))
    np.save(os.path.join(OUT, f"{e}_tf.npy"), T[:, k].astype(np.float32))
    np.save(os.path.join(OUT, f"{e}_rfp.npy"), q.rfp.values[k].astype(np.float32))
    names = feat.id[isgene].tolist()
    names = ["ATMGatp4" if n == "atp4" or f == "atp4" else n for n, f in zip(names, feat.name[isgene])]
    pd.Series(names).to_csv(os.path.join(OUT, f"{e}_genes.txt"), index=False, header=False)
    pd.Series(tf_names).to_csv(os.path.join(OUT, f"{e}_tfnames.txt"), index=False, header=False)
    q.loc[k, ["cell", "library"]].to_csv(os.path.join(OUT, f"{e}_cells.tsv"), sep="\t", index=False)
    s = dict(cells_cellranger=len(q), in_paper_matrix=int(q.in_paper_matrix.sum()), paper_matrix_cells=int(len(paper)),
             low_umi=int(q.low_umi.sum()), low_genes=int(q.low_genes.sum()), high_mito=int(q.high_mito.sum()),
             high_plastid=int(q.high_plastid.sum()), doublets=int(q.doublet.sum()), pass_qc=int(k.sum()),
             umi_median_pass=float(q.umi[k].median()), genes_median_pass=float(q.genes[k].median()),
             frac_no_TF_pass=float((q.tf_reads[k] == 0).mean()), median_nTF_pass=float(q.n_tf[k].median()),
             RFP_detected_pass=float((q.rfp[k] > 0).mean()), GFP_detected=float((q.gfp > 0).mean()), tfs=tf_names)
    summary[e] = s
    print(e, {a: b for a, b in s.items() if a != "tfs"}, flush=True)
json.dump(summary, open(os.path.join(QC, "qc_summary.json"), "w"), indent=1)
