#!/usr/bin/env python3
"""Step 66: are the scTARGET-seq TF effects and TF x TF interactions (Fig 5 network, 461 interactions)
visible with ENDOGENOUS TF expression in an independent atlas (Shahan et al. 2022, Dev Cell; 110,427 WT
root cells, 16 samples; no TF over-expression)?

Model (the scTARGET logic, adapted to endogenous TFs), fitted for the 11,980 scTARGET genes:
  y    = log1p(CP10k) of the gene (raw UMIs, spliced + unspliced; depth = all genes)
  a, b = endogenous detection of TFa / TFb in the cell (>= 1 raw UMI at the TF locus)
  interaction: y ~ a + b + a:b + stratum + covariates
  main effect: y ~ a + stratum + covariates
  stratum    = sample x cell type x developmental time group (Shahan annotation; celltype.anno x
               consensus.time.group), absorbed by within-stratum demeaning (Frisch-Waugh-Lovell)
  covariates = log depth, (log depth)^2, log genes detected, % mitochondrial reads
  HC1 robust standard errors.  Sensitivity: TF dose = log1p(UMI) instead of detection.
Null: for each network edge (TFa, TFb, target), the same pair tested on control genes of the scTARGET
universe that are not targets (any tier) of TFa or TFb in scTARGET, matched on atlas expression (20 bins).
Concordance = sign(atlas estimate) == sign(scTARGET meta-analysis estimate).
"""
import os, sys, json
import numpy as np, pandas as pd
from scipy import sparse, stats
sys.path.insert(0, os.path.dirname(__file__))
from common import RES, TF_AGI, bh

ATLAS = os.path.expanduser("~/scdata/Shahan2022")
OUT = os.path.join(RES, "atlas"); os.makedirs(OUT, exist_ok=True)
AGI = dict(TF_AGI); AGI["HHO3"] = "AT1G25550"      # endogenous HHO3 (NIGT1.1) locus
MIN_STATE = 20                                     # testable pair: >= 20 cells in each of A-only, B-only, both
rng = np.random.default_rng(66)

# ------------------------------------------------------------------ atlas counts and annotation
D = os.path.join(ATLAS, "derived")
ng, nc = [int(v) for v in open(os.path.join(D, "dim.txt")).read().split()]
M = sparse.csc_matrix((np.fromfile(os.path.join(D, "counts_x.int32"), np.int32).astype(np.float32),
                       np.fromfile(os.path.join(D, "counts_i.int32"), np.int32),
                       np.fromfile(os.path.join(D, "counts_p.int32"), np.int32)), shape=(ng, nc))
agenes = pd.Index(open(os.path.join(D, "genes.txt")).read().split())
acells = open(os.path.join(D, "cells.txt")).read().split()
meta = pd.read_csv(os.path.join(ATLAS, "metadata", "Shahan2022_atlas_cell_metadata.tsv.gz"), sep="\t", low_memory=False)
meta = meta.set_index("cell_barcode_id")
# count-matrix cell names carry merge suffixes (AAAC..._1_1_1_1_1_1_1_1) while Data S3 uses AAAC..._1:
# the two tables are in the same order; check barcode + spliced/unspliced UMI totals for every cell.
cc = pd.read_csv(os.path.join(D, "cells_counts.tsv"), sep="\t")
assert len(cc) == len(meta) == M.shape[1]
same = ((cc.cell.str[:16].values == meta.index.str[:16].values)
        & (cc.nCount_spliced.values == meta.nCount_spliced_RNA.astype(int).values)
        & (cc.nCount_unspliced.values == meta.nCount_unspliced_RNA.astype(int).values))
assert same.all(), f"{(~same).sum()} cells do not match"
cid = meta.index
n = M.shape[1]
tot = np.asarray(M.sum(0)).ravel().astype(np.float64)
ngen = np.diff(M.indptr).astype(np.float64)
print(f"atlas: {n} cells x {ng} genes; depth vs annotation r = {np.corrcoef(tot, meta.nCount_RNA.astype(float))[0, 1]:.4f}")

strat = (meta["orig.ident"] + "|" + meta["celltype.anno"] + "|" + meta["consensus.time.group"]).values
codes, levels = pd.factorize(strat); K = len(levels)
cnt = np.bincount(codes, minlength=K).astype(np.float64)
S = sparse.csr_matrix((np.ones(n), (np.arange(n), codes)), shape=(n, K))


def demean(X):
    X = np.asarray(X, np.float64)
    return X - ((S.T @ X) / cnt[:, None])[codes]


C = demean(np.column_stack([np.log(tot), np.log(tot) ** 2, np.log(ngen), meta["percent.mt"].astype(float).values]))
Cp = np.linalg.pinv(C)


def tilde(X):
    X = demean(X)
    return X - C @ (Cp @ X)


DF = n - K - C.shape[1]

# ------------------------------------------------------------------ scTARGET results
SRC = os.path.join(RES, "final")
F = np.load(os.path.join(SRC, "full_stats.npz"))
U = pd.Index(F["genes"]); tfs = list(F["tfs"]); pairs_all = list(F["pairs"])
for t in tfs:
    if AGI[t] not in agenes: print("TF locus missing in atlas:", t, AGI[t])


def meta_fe(kind):
    b1, s1, b2, s2 = (F[f"{e}_{kind}_{v}"] for e in ("Exp1", "Exp2") for v in ("beta", "se"))
    w1, w2 = 1 / s1 ** 2, 1 / s2 ** 2
    b = (w1 * b1 + w2 * b2) / (w1 + w2)
    return b, b * np.sqrt(w1 + w2)


sc_main_b, sc_main_z = meta_fe("main"); sc_int_b, sc_int_z = meta_fe("inter")
NET = pd.read_csv(os.path.join(RES, "paper", "fig5d_candidates.tsv"), sep="\t")
IE = pd.read_csv(os.path.join(SRC, "inter_effects.tsv.gz"), sep="\t")
ME = pd.read_csv(os.path.join(SRC, "main_effects.tsv.gz"), sep="\t")
NET = NET.merge(IE[["pair", "target", "beta_meta", "q_meta"]].rename(columns={"q_meta": "q_meta_ie"}), on=["pair", "target"], how="left")
NET["a"] = NET.pair.str.split("-").str[0]; NET["b"] = NET.pair.str.split("-").str[1]
MAIN = ME[ME.tier.isin(["strict", "medium"])].rename(columns={"regulator": "tf"}).copy()
pairs = sorted(NET.pair.unique())
print(f"network: {len(NET)} interactions, {len(pairs)} pairs, {NET.target.nunique()} targets; main effects: {len(MAIN)}")

# genes of the scTARGET universe present in the atlas
G = [g for g in U if g in agenes and g not in set(AGI.values())]
gi_atlas = agenes.get_indexer(G)
Msub = M[gi_atlas, :].tocsr()
print(f"universe genes in atlas: {len(G)} / {len(U)}")

# ------------------------------------------------------------------ TF states
def tf_row(t):
    return np.asarray(M[agenes.get_loc(AGI[t]), :].todense()).ravel()


TFc = pd.DataFrame({t: tf_row(t) for t in tfs if AGI[t] in agenes})
DET = (TFc >= 1).astype(np.float64)
DOSE = np.log1p(TFc)
det_ct = DET.groupby(meta["celltype.anno"].values).mean()
det_ct.loc["All cells"] = DET.mean()
det_ct.to_csv(os.path.join(OUT, "tf_detection_by_celltype.tsv"), sep="\t")
det_tz = DET.groupby((meta["celltype.anno"] + "|" + meta["time.anno"]).values).mean()
det_tz.to_csv(os.path.join(OUT, "tf_detection_by_celltype_zone.tsv"), sep="\t")
print("TF detection (% cells):", (100 * DET.mean()).round(2).sort_values().to_dict())

# design matrices (residualised): main = 1 column per TF; interaction = 3 columns per pair
Xm = tilde(DET.values); Xm_d = tilde(DOSE.values)
pc = []
for p in pairs:
    a, b = p.split("-")
    for V in (DET, DOSE):
        pc.append(np.column_stack([V[a], V[b], V[a] * V[b]]))
cnt4 = []
for p in pairs:
    a, b = p.split("-")
    A_, B_ = DET[a].values > 0, DET[b].values > 0
    both_strata = len(set(codes[A_ & B_]))
    cnt4.append(dict(pair=p, n_neither=int((~A_ & ~B_).sum()), n_A=int((A_ & ~B_).sum()), n_B=int((~A_ & B_).sum()),
                     n_both=int((A_ & B_).sum()), strata_with_both=both_strata,
                     expected_both=float(((S.T @ A_.astype(float)) * (S.T @ B_.astype(float)) / cnt).sum())))
CNT4 = pd.DataFrame(cnt4)
CNT4["testable"] = (CNT4[["n_A", "n_B", "n_both"]] >= MIN_STATE).all(1)
CNT4["codetect_ratio"] = CNT4.n_both / CNT4.expected_both
CNT4.to_csv(os.path.join(OUT, "pair_state_counts.tsv"), sep="\t", index=False)
print(f"testable pairs: {CNT4.testable.sum()} / {len(CNT4)}")


# stacked operators for all pairs x {detection, dose}: beta = A @ y; robust variance by expansion
# var_k = sum_i A_ki^2 e_i^2, e = y - X beta  ->  T1 - 2 sum_l b_l T2_kl + sum_lm b_l b_m T3_klm  (all GEMMs)
nP = len(pairs)
Xall = np.empty((nP, 2, n, 3)); Aall = np.empty((nP, 2, 3, n))
for j in range(nP):
    for v in range(2):
        X3 = tilde(pc[2 * j + v]); Xall[j, v] = X3
        Aall[j, v] = np.linalg.solve(X3.T @ X3, X3.T)
del pc
W = Aall ** 2                                                        # nP,2,3(k),n
T3 = np.einsum("pvkn,pvnl,pvnm->pvklm", W, Xall, Xall)
WX = np.einsum("pvkn,pvnl->pvkln", W, Xall)                          # nP,2,3,3,n
A_flat = Aall.reshape(-1, n); W_flat = W.reshape(-1, n); WX_flat = WX.reshape(-1, n)
del WX, Xall

# ------------------------------------------------------------------ fit all genes in chunks
nG = len(G); CH = 300
res = {k: np.full((nP, nG), np.nan) for k in ("int_b", "int_se", "int_b_dose", "int_se_dose", "bA", "bB", "seA", "seB")}
main = {k: np.full((Xm.shape[1], nG), np.nan) for k in ("b", "se", "b_dose", "se_dose")}
mean_expr = np.zeros(nG); det_frac = np.zeros(nG)
sxx = (Xm ** 2).sum(0); sxx_d = (Xm_d ** 2).sum(0)
for s0 in range(0, nG, CH):
    blk = Msub[s0:s0 + CH].toarray().T.astype(np.float64)            # n x g raw counts
    g_ = blk.shape[1]; sl = slice(s0, s0 + g_)
    det_frac[sl] = (blk > 0).mean(0)
    Y = np.log1p(blk / tot[:, None] * 1e4); del blk
    mean_expr[sl] = Y.mean(0)
    Yt = tilde(Y); Y2 = Yt ** 2; del Y
    for X, sx, kb, ks in ((Xm, sxx, "b", "se"), (Xm_d, sxx_d, "b_dose", "se_dose")):
        b = (X.T @ Yt) / sx[:, None]
        v = (X.T ** 2 @ Y2) - 2 * b * ((X.T ** 3) @ Yt) + b ** 2 * (X ** 4).sum(0)[:, None]
        main[kb][:, sl] = b; main[ks][:, sl] = np.sqrt(np.maximum(v, 0) / sx[:, None] ** 2 * n / DF)
    B = (A_flat @ Yt).reshape(nP, 2, 3, g_)
    T1 = (W_flat @ Y2).reshape(nP, 2, 3, g_)
    T2 = (WX_flat @ Yt).reshape(nP, 2, 3, 3, g_)
    V = T1 - 2 * np.einsum("pvlg,pvklg->pvkg", B, T2) + np.einsum("pvlg,pvmg,pvklm->pvkg", B, B, T3)
    SE = np.sqrt(np.maximum(V, 0) * n / (DF - 3))
    res["int_b"][:, sl] = B[:, 0, 2]; res["int_se"][:, sl] = SE[:, 0, 2]
    res["int_b_dose"][:, sl] = B[:, 1, 2]; res["int_se_dose"][:, sl] = SE[:, 1, 2]
    res["bA"][:, sl] = B[:, 0, 0]; res["bB"][:, sl] = B[:, 0, 1]; res["seA"][:, sl] = SE[:, 0, 0]; res["seB"][:, sl] = SE[:, 0, 1]
    print(f"  genes {s0 + g_}/{nG}", flush=True)
np.savez_compressed(os.path.join(OUT, "atlas_full_stats.npz"), genes=np.array(G), pairs=np.array(pairs),
                    tfs=np.array(DET.columns), mean_expr=mean_expr, det_frac=det_frac,
                    **{f"int_{k}": v for k, v in res.items()}, **{f"main_{k}": v for k, v in main.items()})

# ------------------------------------------------------------------ edge-level tests vs matched control genes
gpos = pd.Index(G)
ebin = pd.qcut(pd.Series(mean_expr).rank(method="first"), 20, labels=False).values
sc_gidx = U.get_indexer(G)
tf_targets = {t: set(ME.loc[ME.regulator == t, "target"]) | set(IE.loc[IE.pair.str.split("-").str[0].eq(t) | IE.pair.str.split("-").str[1].eq(t), "target"]) for t in tfs}


def edge_table(rows, est, se, groups, eligible):
    """rows: DataFrame with 'gene'; est/se: arrays (n_groups x nG); groups: row index into est; eligible: per-row gene mask."""
    out = []
    for r, g_row, ok in zip(rows.itertuples(), groups, eligible):
        gi = gpos.get_loc(r.gene)
        z = est[g_row, gi] / se[g_row, gi]
        ctrl = np.where(ok & (ebin == ebin[gi]))[0]
        zc = est[g_row, ctrl] / se[g_row, ctrl]
        zc = zc[np.isfinite(zc)]
        s = r.ref_sign
        out.append(dict(est=est[g_row, gi], se=se[g_row, gi], z=z, p=2 * stats.norm.sf(abs(z)),
                        concordant=bool(np.sign(z) == s), z_signed=z * s,
                        p_emp=(1 + (zc * s >= z * s).sum()) / (1 + len(zc)), n_ctrl=len(zc),
                        ctrl_z_mean=zc.mean(), ctrl_z_sd=zc.std(), ctrl_conc=(np.sign(zc) == s).mean()))
    return pd.concat([rows.reset_index(drop=True), pd.DataFrame(out)], axis=1)


def null_draws(T, est, se, groups, eligible, R=2000):
    """Concordance and n(p<0.05, same sign) when each edge's target is replaced by one random matched control."""
    conc = np.zeros(R); nsig = np.zeros(R)
    picks = []
    for r, g_row, ok in zip(T.itertuples(), groups, eligible):
        gi = gpos.get_loc(r.gene)
        ctrl = np.where(ok & (ebin == ebin[gi]) & np.isfinite(est[g_row]) & (se[g_row] > 0))[0]
        zc = est[g_row, ctrl] / se[g_row, ctrl]
        picks.append(zc[rng.integers(0, len(zc), R)] * r.ref_sign)
    Zs = np.array(picks)                                               # edges x R (signed)
    conc = (Zs > 0).mean(0); nsig = ((Zs > 0) & (2 * stats.norm.sf(np.abs(Zs)) < 0.05)).sum(0)
    return conc, nsig, Zs


summary = {}
# interactions (network)
pair_row = {p: i for i, p in enumerate(pairs)}
T = NET.rename(columns={"target": "gene"}).copy()
T = T[T.gene.isin(gpos)].copy()
T["ref_sign"] = np.sign(T.beta_meta.fillna(T.nonadd))
T = T.merge(CNT4[["pair", "n_A", "n_B", "n_both", "testable", "codetect_ratio"]], on="pair", how="left")
Garr = np.array(G)
pair_mask = {p: ~np.isin(Garr, list(tf_targets[p.split("-")[0]] | tf_targets[p.split("-")[1]])) for p in pairs}
tf_mask = {t: ~np.isin(Garr, list(tf_targets[t])) for t in tfs}
elig = [pair_mask[p] for p in T.pair]
grp = [pair_row[p] for p in T.pair]
TI = edge_table(T, res["int_b"], res["int_se"], grp, elig)
TI_d = edge_table(T, res["int_b_dose"], res["int_se_dose"], grp, elig)
TI["z_dose"] = TI_d.z.values; TI["concordant_dose"] = TI_d.concordant.values
TI["bA"] = [res["bA"][g, gpos.get_loc(x)] for g, x in zip(grp, T.gene)]
TI["bB"] = [res["bB"][g, gpos.get_loc(x)] for g, x in zip(grp, T.gene)]
TI["atlas_mean_expr"] = [mean_expr[gpos.get_loc(x)] for x in T.gene]
TI["atlas_det_frac"] = [det_frac[gpos.get_loc(x)] for x in T.gene]
tt = TI.testable.values
TI["q_emp"] = np.nan; TI.loc[tt, "q_emp"] = bh(TI.loc[tt, "p_emp"].values)
TI["sc_z"] = [sc_int_z[pairs_all.index(p), U.get_loc(g)] for p, g in zip(TI.pair, TI.gene)]
TI.to_csv(os.path.join(OUT, "network_edges_atlas.tsv"), sep="\t", index=False)
conc_null, nsig_null, _ = null_draws(T[tt], res["int_b"], res["int_se"], [g for g, k in zip(grp, tt) if k], [e for e, k in zip(elig, tt) if k])
obs_c = TI.loc[tt, "concordant"].mean(); obs_s = int(((TI.z_signed > 0) & (TI.p < 0.05))[tt].sum())
summary["interactions"] = dict(
    n_network=len(NET), n_in_atlas=len(TI), n_testable=int(tt.sum()), n_pairs_testable=int(CNT4.testable.sum()),
    concordance=obs_c, concordance_null_mean=conc_null.mean(), concordance_null_95=float(np.quantile(conc_null, 0.95)),
    p_concordance=float((1 + (conc_null >= obs_c).sum()) / (1 + len(conc_null))),
    n_same_sign_p05=obs_s, n_same_sign_p05_null_mean=float(nsig_null.mean()), n_same_sign_p05_null_95=float(np.quantile(nsig_null, 0.95)),
    p_n_same_sign_p05=float((1 + (nsig_null >= obs_s).sum()) / (1 + len(nsig_null))),
    n_q_emp_lt_0_1=int((TI.q_emp < 0.1).sum()),
    spearman_sc_vs_atlas=float(stats.spearmanr(TI.loc[tt, "sc_z"], TI.loc[tt, "z"]).correlation),
    concordance_dose=float(TI.loc[tt, "concordant_dose"].mean()),
    by_shape={s: dict(n=int(((TI["shape"] == s) & tt).sum()), concordance=float(TI.loc[(TI["shape"] == s) & tt, "concordant"].mean()))
              for s in TI["shape"].unique()})
np.save(os.path.join(OUT, "null_interaction_concordance.npy"), conc_null)

# main effects (positive control): scTARGET TF -> target, endogenous TF detection
tf_col = {t: i for i, t in enumerate(DET.columns)}
Tm = MAIN.rename(columns={"target": "gene"})
Tm = Tm[Tm.gene.isin(gpos) & Tm.tf.isin(tf_col)].copy()
Tm["ref_sign"] = np.sign(Tm.beta_meta)
eligm = [tf_mask[t] for t in Tm.tf]
grpm = [tf_col[t] for t in Tm.tf]
TM = edge_table(Tm, main["b"], main["se"], grpm, eligm)
TM["sc_z"] = [sc_main_z[tfs.index(t), U.get_loc(g)] for t, g in zip(TM.tf, TM.gene)]
TM["q_emp"] = bh(TM.p_emp.values)
TM.to_csv(os.path.join(OUT, "main_effects_atlas.tsv"), sep="\t", index=False)
conc_null_m, nsig_null_m, _ = null_draws(Tm, main["b"], main["se"], grpm, eligm)
obs_cm = TM.concordant.mean(); obs_sm = int(((TM.z_signed > 0) & (TM.p < 0.05)).sum())
summary["main_effects"] = dict(
    n=len(TM), concordance=obs_cm, concordance_null_mean=conc_null_m.mean(), concordance_null_95=float(np.quantile(conc_null_m, 0.95)),
    p_concordance=float((1 + (conc_null_m >= obs_cm).sum()) / (1 + len(conc_null_m))),
    n_same_sign_p05=obs_sm, n_same_sign_p05_null_mean=float(nsig_null_m.mean()),
    p_n_same_sign_p05=float((1 + (nsig_null_m >= obs_sm).sum()) / (1 + len(nsig_null_m))),
    n_q_emp_lt_0_1=int((TM.q_emp < 0.1).sum()),
    spearman_sc_vs_atlas=float(stats.spearmanr(TM.sc_z, TM.z).correlation),
    by_tf=TM.groupby("tf").concordant.agg(["size", "mean"]).round(3).to_dict("index"))
np.save(os.path.join(OUT, "null_main_concordance.npy"), conc_null_m)

summary["atlas"] = dict(cells=int(n), strata=int(K), samples=int(meta["orig.ident"].nunique()),
                        tf_detection_pct=(100 * DET.mean()).round(3).to_dict(), min_state=MIN_STATE)
json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=1, default=float)
print(json.dumps(summary, indent=1, default=float))
