#!/usr/bin/env python3
"""Step 69: the paper's network connections (Fig 5B) in the Shahan 2022 atlas, (i) cell type by cell type and
(ii) with less (or no) cell-type information in the adjustment.

Same per-connection model as step 67 (four endogenous states, gate = dAB - dA - dB, HC1), fitted under:
  global variants (all 110,427 cells)
    full          stratum = sample x cell type x developmental time group   (= step 67)
    no_celltype   stratum = sample x developmental time group (cell type removed)
    sample_only   stratum = sample (cell type and developmental time removed)
  per cell type   cells of one annotated cell type, stratum = sample x developmental time group
A connection is testable in a cell subset when each of TFa-only, TFb-only and both has >= 20 cells.
Recovered = gate significant (BH q < 0.05 within the variant; for cell types, over all connection x cell-type
tests jointly) with the same sign as in scTARGET-seq.
Chance expectation: the scTARGET labels (gate sign and type) are permuted among the targets of the same TF pair
(atlas fits unchanged), 2,000 permutations.
"""
import os, sys, json
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from common import RES, bh
from atlas_common import load_atlas, make_tilde, AGI, STRATA_FULL

OUT = os.path.join(RES, "atlas"); os.makedirs(OUT, exist_ok=True)
MIN_STATE = 20
rng = np.random.default_rng(69)
NET = pd.read_csv(os.path.join(RES, "paper", "fig5d_candidates.tsv"), sep="\t")
NET["gate"] = NET["shape"].replace({"cancellation": "competition"})
A = load_atlas(); meta, tot, ngen = A["meta"], A["tot"], A["ngen"]
tfs = sorted(set(NET.pair.str.split("-").str[0]) | set(NET.pair.str.split("-").str[1]))
need = [g for g in sorted(set(NET.target) | {AGI[t] for t in tfs}) if g in A["genes"]]
R = A["M"][A["genes"].get_indexer(need), :].tocsr()
row = {g: np.asarray(R[i].todense()).ravel().astype(np.float64) for i, g in enumerate(need)}
del R, A
NET = NET[NET.target.isin(row)].reset_index(drop=True)
det = {t: row[AGI[t]] >= 1 for t in tfs}
logy = {g: np.log1p(row[g] / tot * 1e4) for g in set(NET.target)}


def fit(mask, keys, label):
    tilde, DF = make_tilde(meta, tot, ngen, mask, keys)
    n = int(mask.sum()); ycache = {}; out = []
    for p, sub in NET.groupby("pair"):
        a, b = p.split("-")
        sa, sb = det[a][mask], det[b][mask]
        S3 = np.column_stack([sa & ~sb, ~sa & sb, sa & sb]).astype(np.float64)
        ns = S3.sum(0)
        if (ns < MIN_STATE).any():
            continue
        X = tilde(S3); XtXi = np.linalg.inv(X.T @ X)
        for r in sub.itertuples():
            if r.target not in ycache:
                ycache[r.target] = tilde(logy[r.target][mask][:, None])[:, 0]
            y = ycache[r.target]
            if not np.any(y):
                continue
            beta = XtXi @ (X.T @ y); e = y - X @ beta
            V = XtXi @ ((X * (e ** 2)[:, None]).T @ X) @ XtXi * n / (DF - 3)
            c = np.array([-1.0, -1.0, 1.0]); na = c @ beta; se = np.sqrt(c @ V @ c); sd = np.sqrt(np.diag(V))
            out.append(dict(variant=label, pair=p, target=r.target, gate_paper=r.gate, paper_nonadd=r.nonadd,
                            paper_dA=r.dA, paper_dB=r.dB, paper_dAB=r.dAB, n_cells=n, n_A=int(ns[0]), n_B=int(ns[1]), n_both=int(ns[2]),
                            dA=beta[0], dB=beta[1], dAB=beta[2], se_dA=sd[0], se_dB=sd[1], se_dAB=sd[2], p_dA=2 * stats.norm.sf(abs(beta[0] / sd[0])),
                            p_dB=2 * stats.norm.sf(abs(beta[1] / sd[1])), nonadd=na, se_nonadd=se, z_nonadd=na / se,
                            p_nonadd=2 * stats.norm.sf(abs(na / se))))
    print(f"  {label}: {n} cells, {len(out)} connections testable", flush=True)
    return pd.DataFrame(out)


def classify(T):
    sA, sB = T.p_dA < 0.05, T.p_dB < 0.05
    s1 = np.sign(np.where(sA, T.dA, 0) + np.where(sB, T.dB, 0))
    return np.select([~sA & ~sB, np.sign(T.dAB) == -s1, np.sign(T.nonadd) == s1], ["emergent", "reversal", "synergy"], "competition")


ALL = np.ones(len(meta), bool)
variants = {"full": STRATA_FULL, "no_celltype": ["orig.ident", "consensus.time.group"], "sample_only": ["orig.ident"]}
res = [fit(ALL, k, v) for v, k in variants.items()]
cts = meta["celltype.anno"].value_counts().index.tolist()
for ct in cts:
    res.append(fit((meta["celltype.anno"] == ct).values, ["orig.ident", "consensus.time.group"], ct))
T = pd.concat(res, ignore_index=True)
T["is_celltype"] = ~T.variant.isin(list(variants))
T["q_nonadd"] = np.nan
for v in variants:
    m = T.variant == v; T.loc[m, "q_nonadd"] = bh(T.loc[m, "p_nonadd"].values)
m = T.is_celltype; T.loc[m, "q_nonadd"] = bh(T.loc[m, "p_nonadd"].values)
T["same_sign"] = np.sign(T.nonadd) == np.sign(T.paper_nonadd)
T["status"] = np.where(T.q_nonadd >= 0.05, "not significant", np.where(T.same_sign, "recovered", "opposite"))
T["gate_atlas"] = np.where(T.q_nonadd < 0.05, classify(T), "none")
T["same_gate"] = T.gate_atlas == T.gate_paper
T.to_csv(os.path.join(OUT, "connections_by_variant_and_celltype.tsv.gz"), sep="\t", index=False)

# ---------------------------------------------------------------- chance expectation by permutation of scTARGET labels
conn = NET[["pair", "target", "nonadd", "gate"]].rename(columns={"nonadd": "paper_nonadd", "gate": "gate_paper"})
groups = [np.where(conn.pair.values == p)[0] for p in conn.pair.unique()]
key = conn.pair + "|" + conn.target
T["key"] = T.pair + "|" + T.target
pos = pd.Series(np.arange(len(conn)), index=key.values)
T["ci"] = pos.loc[T.key].values
sig = T.q_nonadd.values < 0.05; sgn = np.sign(T.nonadd.values); ga = T.gate_atlas.values; ci = T.ci.values


def stats_for(psign, pgate):
    rec = sig & (sgn == psign[ci]); recg = rec & (ga == pgate[ci]); opp = sig & (sgn != psign[ci])
    out = {}
    for v in list(variants):
        m = (T.variant == v).values
        out[v] = (rec[m].sum(), recg[m].sum(), opp[m].sum())
    m = T.is_celltype.values
    out["celltypes"] = (len(set(ci[m & rec])), len(set(ci[m & recg])), len(set(ci[m & opp])), rec[m].sum(), opp[m].sum())
    return out


obs = stats_for(np.sign(conn.paper_nonadd.values), conn.gate_paper.values)
null = []
for _ in range(2000):
    i = np.arange(len(conn))
    for g in groups:
        i[g] = rng.permutation(g)
    null.append(stats_for(np.sign(conn.paper_nonadd.values)[i], conn.gate_paper.values[i]))

summ = {}
names3 = ["recovered", "recovered_same_gate", "opposite"]
for v in variants:
    d = T[T.variant == v]
    s = {"n_testable": len(d), "status": d.status.value_counts().to_dict(),
         "rho_nonadd": float(stats.spearmanr(d.nonadd, d.paper_nonadd).correlation),
         "rho_dAB": float(stats.spearmanr(d.dAB, d.paper_dAB).correlation)}
    for j, nm in enumerate(names3):
        nv = np.array([x[v][j] for x in null])
        s[nm] = dict(obs=int(obs[v][j]), chance_mean=float(nv.mean()), chance_95=[float(np.percentile(nv, 2.5)), float(np.percentile(nv, 97.5))],
                     p_more_than_chance=float((1 + (nv >= obs[v][j]).sum()) / (1 + len(nv))))
    summ[v] = s
names5 = ["connections_recovered_in_>=1_celltype", "connections_recovered_same_gate_in_>=1_celltype",
          "connections_opposite_in_>=1_celltype", "tests_recovered", "tests_opposite"]
c = T[T.is_celltype]
summ["per_celltype"] = {"tests": len(c), "connections_testable_in_>=1_celltype": int(c.ci.nunique()),
                        "by_celltype": c.groupby("variant").apply(lambda d: {"testable": len(d), "recovered": int((d.status == "recovered").sum()),
                                                                             "opposite": int((d.status == "opposite").sum())}).to_dict()}
for j, nm in enumerate(names5):
    nv = np.array([x["celltypes"][j] for x in null])
    summ["per_celltype"][nm] = dict(obs=int(obs["celltypes"][j]), chance_mean=float(nv.mean()),
                                    chance_95=[float(np.percentile(nv, 2.5)), float(np.percentile(nv, 97.5))],
                                    p_more_than_chance=float((1 + (nv >= obs["celltypes"][j]).sum()) / (1 + len(nv))))
# per cell type (and per global variant): gates recovered / opposite and the combined effect (both TFs vs neither)
# against the same within-pair permutation of the scTARGET labels
pc = {}
for v in list(variants) + list(c.variant.unique()):
    d = T[T.variant == v].reset_index(drop=True)
    gs = [np.where(d.pair.values == p)[0] for p in d.pair.unique()]
    sg = d.q_nonadd.values < 0.05; sn = np.sign(d.nonadd.values)
    def one(i):
        ps = np.sign(d.paper_nonadd.values[i])
        return ((sg & (sn == ps)).sum(), (sg & (sn != ps)).sum(), stats.spearmanr(d.dAB, d.paper_dAB.values[i]).correlation,
                stats.spearmanr(d.nonadd, d.paper_nonadd.values[i]).correlation)
    o = one(np.arange(len(d)))
    nv = []
    for _ in range(1000):
        i = np.arange(len(d))
        for g in gs:
            i[g] = rng.permutation(g)
        nv.append(one(i))
    nv = np.array(nv)
    pc[v] = dict(testable=len(d), recovered=int(o[0]), recovered_chance=float(nv[:, 0].mean()),
                 p_recovered=float((1 + (nv[:, 0] >= o[0]).sum()) / 1001),
                 opposite=int(o[1]), opposite_chance=float(nv[:, 1].mean()),
                 rho_both=float(o[2]), rho_both_chance=float(nv[:, 2].mean()), p_rho_both=float((1 + (nv[:, 2] >= o[2]).sum()) / 1001),
                 rho_gate=float(o[3]), rho_gate_chance=float(nv[:, 3].mean()), p_rho_gate=float((1 + (nv[:, 3] >= o[3]).sum()) / 1001),
                 recovered_same_gate=int(((d.status == "recovered") & d.same_gate).sum()))
    print(v, pc[v], flush=True)
summ["per_celltype"]["recovered_vs_chance"] = pc
pd.DataFrame([dict(variant=k, **v) for k, v in pc.items()]).to_csv(os.path.join(OUT, "celltype_and_variant_summary.tsv"), sep="\t", index=False)
json.dump(summ, open(os.path.join(OUT, "celltype_strata_summary.json"), "w"), indent=1, default=float)
print(json.dumps(summ, indent=1, default=float))
