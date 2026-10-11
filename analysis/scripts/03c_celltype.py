#!/usr/bin/env python3
"""Step 3c: marker-based root cell identity, TF distribution across cell
types, and the effect of adjusting for cell identity on reproducibility.

Cell-type scores = mean z-score (log1p CP10k) of literature marker genes.
Each cell gets the type with the highest score (if > 0.5 SD), else
'unassigned'.  Scores are independent of the TF design by construction
(markers, not data-driven PCs), so adjusting for them cannot absorb a TF's
genome-wide program unless the TF regulates the markers themselves.
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from common import *

MARKERS = {
    "root_hair": ["AT5G49270", "AT1G12560", "AT1G62980", "AT3G62680"],
    "non_hair_epidermis": ["AT1G79840", "AT5G14750"],
    "cortex": ["AT1G62510", "AT1G62500"],
    "endodermis": ["AT3G54220", "AT5G57620", "AT2G36100"],
    "xylem": ["AT4G35350", "AT1G20850", "AT5G17420", "AT4G18780", "AT5G44030", "AT1G71930"],
    "phloem": ["AT1G22710", "AT1G79430", "AT3G01680", "AT4G19840"],
    "root_cap_columella": ["AT1G79580", "AT1G33280", "AT4G10350", "AT1G26870", "AT1G70940"],
    "meristem_dividing": ["AT4G37490", "AT2G28740", "AT1G08560", "AT3G20840", "AT1G51190"],
    "stele_SHR_QC": ["AT4G37650", "AT3G11260"],
}

ds = {e: load(e) for e in EXPS}
genes, tfs = analysis_sets(ds)
rows, ct_tables = [], {}
for e in EXPS:
    d = ds[e]
    G = d["G"]; tot = G.sum(0).A1
    gi = pd.Index(d["genes"])
    S = {}
    for ct, mk in MARKERS.items():
        idx = [gi.get_loc(g) for g in mk if g in gi]
        if not idx:
            continue
        L = np.log1p(G[idx, :].toarray().T / tot[:, None] * 1e4)
        Zs = (L - L.mean(0)) / (L.std(0) + 1e-9)
        S[ct] = Zs.mean(1)
        rows.append(dict(exp=e, cell_type=ct, markers_found=len(idx), markers=",".join(np.array(mk)[[g in gi for g in mk]])))
    S = pd.DataFrame(S)
    call = np.where(S.max(1) > 0.5, S.idxmax(1), "unassigned")
    S["cell_type"] = call
    S["cell"] = d["cells"]["cell"].values
    S.to_csv(os.path.join(RES, f"{e}_celltype_scores.tsv"), sep="\t", index=False)
    vc = pd.Series(call).value_counts()
    print(f"\n{e} cell-type calls:", vc.to_dict())
    # TF distribution across cell types: % of cells of each type with TF (>=1)
    P = (d["T"][tfs].values >= 1)
    tab = pd.DataFrame(P, columns=tfs).groupby(call).mean().T * 100
    tab.loc["n_cells"] = vc.reindex(tab.columns).values
    tab.loc["mean_n_TFs"] = pd.Series(P.sum(1)).groupby(call).mean().reindex(tab.columns).values
    ct_tables[e] = tab
    tab.round(1).to_csv(os.path.join(RES, f"{e}_TF_by_celltype.tsv"), sep="\t")
    print(tab.round(1).to_string())
    # chi-square per TF (presence x cell type), adjusted for nothing
    chi = {}
    for j, tf in enumerate(tfs):
        ct = pd.crosstab(call, P[:, j])
        chi[tf] = stats.chi2_contingency(ct)[1]
    print("TF x cell-type association p-values (chi2):", {k: f"{v:.1e}" for k, v in chi.items()})
    # does the TF design predict cell-type scores beyond load?  (R2 increase)
    X = encode(d["T"], tfs, "ge1"); Z = covariates(d, tfs).values
    D0 = np.column_stack([np.ones(len(X)), Z]); D1 = np.column_stack([D0, X])
    for ct in MARKERS:
        if ct not in S: continue
        y = S[ct].values
        r0 = y - D0 @ np.linalg.lstsq(D0, y, rcond=None)[0]
        r1 = y - D1 @ np.linalg.lstsq(D1, y, rcond=None)[0]
        n, k = len(y), X.shape[1]
        F = ((r0 ** 2).sum() - (r1 ** 2).sum()) / k / ((r1 ** 2).sum() / (n - D1.shape[1]))
        print(f"  {ct:20s} partial R2 of TF identity = {1 - (r1**2).sum()/(r0**2).sum():.4f}  F-test p = {stats.f.sf(F, k, n - D1.shape[1]):.1e}")

pd.DataFrame(rows).to_csv(os.path.join(RES, "celltype_markers_used.tsv"), sep="\t", index=False)
