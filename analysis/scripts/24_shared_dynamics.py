#!/usr/bin/env python3
"""Step 24: shared ERF5 x CRF4 response *dynamics* across the 4 datasets
(bulk DoubleTARGET exp 1 & 2, single-cell Exp1 & Exp2), beyond sign agreement.

For each gene and dataset: profile = standardized effects (effect / SE) of
ERF5 alone, CRF4 alone and ERF5+CRF4 vs EV/none (single cell: four-state model
adjusted for load, depth, cell type and the other TFs; bulk: replicate means,
pooled SD).  Concordance = mean of the 6 pairwise Pearson correlations of these
profiles between datasets (shape agreement, scale-free), multiplied by the
mean profile strength so flat/noisy genes do not score.  Null: gene labels
permuted independently in every dataset (200x) -> empirical p, BH q.
Concordant genes are then clustered on their concatenated, scaled profiles
into shared dynamic modules."""
import os, sys, re
import numpy as np, pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster
sys.path.insert(0, os.path.dirname(__file__))
from common import *

OUT = os.path.join(RES, "shared_dynamics"); os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(11)
STATES = ["ERF5", "CRF4", "both"]

# ---- bulk
Cn = pd.read_csv("/Users/gabkrouk/Library/CloudStorage/Dropbox/Science/Papers/Emanuel New TARGET paper/Figure Will/counts.txt", sep=" ", index_col=0).drop(columns=["E2_EC.2", "E2_E.1"])
cpm = Cn / Cn.sum() * 1e6
L = np.log2(cpm[(cpm > 1).sum(1) >= 3] + 1)
prof, eff = {}, {}
for e, name in (("E1", "Bulk exp 1"), ("E2", "Bulk exp 2")):
    g = {k: [c for c in L.columns if re.match(e + "_" + k + r"\.", c)] for k in ("N", "E", "C", "EC")}
    m = {k: L[v].mean(1) for k, v in g.items()}
    resid = pd.concat([L[v].sub(m[k], axis=0) for k, v in g.items()], axis=1)
    sd = np.sqrt((resid ** 2).sum(1) / (resid.shape[1] - 4)) + 0.05       # pooled SD, small floor
    E_ = pd.DataFrame({"ERF5": m["E"] - m["N"], "CRF4": m["C"] - m["N"], "both": m["EC"] - m["N"]})
    se = pd.DataFrame({s: sd * np.sqrt(1 / len(g[k]) + 1 / len(g["N"])) for s, k in zip(STATES, ("E", "C", "EC"))})
    eff[name] = E_; prof[name] = E_ / se

# ---- single cell
ds = {e: load(e) for e in EXPS}
genes_sc, tfs = analysis_sets(ds)
CT = {e: pd.read_csv(os.path.join(RES, f"{e}_celltype_scores.tsv"), sep="\t").drop(columns=["cell_type", "cell"]) for e in EXPS}
a, b = tfs.index("ERF5"), tfs.index("CRF4")
for e in EXPS:
    Y = response(ds[e], genes_sc)
    X = encode(ds[e]["T"], tfs, "ge1"); xa, xb = X[:, a], X[:, b]
    S = np.column_stack([xa * (1 - xb), (1 - xa) * xb, xa * xb])
    Z = np.column_stack([covariates(ds[e], tfs, flexible_load=True).values, CT[e].values])
    f = ols(S, np.column_stack([np.delete(X, [a, b], 1), Z]), Y)
    eff[f"Single-cell {e}"] = pd.DataFrame(f["beta"].T, index=genes_sc, columns=STATES)
    prof[f"Single-cell {e}"] = pd.DataFrame(f["t"].T, index=genes_sc, columns=STATES)

names = ["Bulk exp 1", "Bulk exp 2", "Single-cell Exp1", "Single-cell Exp2"]
genes = sorted(set.intersection(*[set(p.index) for p in prof.values()]))
P = np.stack([prof[n].loc[genes, STATES].values for n in names])          # 4 x G x 3
print("genes in all 4 datasets:", len(genes))


def concordance(P):
    Pc = P - P.mean(2, keepdims=True)
    Pn = Pc / (np.linalg.norm(Pc, axis=2, keepdims=True) + 1e-12)
    r = [np.sum(Pn[i] * Pn[j], 1) for i in range(4) for j in range(i + 1, 4)]
    shape = np.mean(r, 0)
    strength = np.sqrt((P ** 2).sum(2)).mean(0)     # mean |t| profile norm
    return shape, shape * np.log1p(strength), np.array(r)


shape, score, R = concordance(P)
null = []
for _ in range(200):
    Pp = np.stack([P[0]] + [P[i][rng.permutation(len(genes))] for i in range(1, 4)])
    null.append(concordance(Pp)[1])
null = np.concatenate(null)
pemp = (1 + (null[None, :] >= score[:, None]).sum(1)) / (1 + null.size) if False else \
       (1 + np.searchsorted(np.sort(null), score, side="left").__rsub__(null.size)) / (1 + null.size)
q = bh(pemp)
res = pd.DataFrame({"gene": genes, "shape_concordance": shape, "score": score, "p_perm": pemp, "q": q})
for k, (i, j) in enumerate([(i, j) for i in range(4) for j in range(i + 1, 4)]):
    res[f"r_{names[i]}__{names[j]}"] = R[k]
for n in names:
    for s in STATES:
        res[f"{n}|{s}"] = eff[n].loc[genes, s].values
res = res.sort_values("score", ascending=False)
for thr in (0.05, 0.1, 0.2):
    print(f"q<{thr}: {(res.q < thr).sum()} genes (expected false ~{thr * (res.q < thr).sum():.0f})")
print("null 95th pct of score:", np.quantile(null, .95).round(3), "| genes above:", (score > np.quantile(null, .95)).sum(), "expected", int(.05 * len(genes)))

# ---- shared dynamic modules among concordant genes
sel = res[res.q < 0.1].copy()
if len(sel) < 30:
    sel = res.head(150).copy()
M = np.stack([prof[n].loc[sel.gene, STATES].values for n in names])     # 4 x g x 3
M = np.concatenate([np.zeros((4, M.shape[1], 1)), M], 2)               # add EV = 0
M = M / (np.abs(M).max(2, keepdims=True) + 1e-12)                       # scale per dataset
F = M.transpose(1, 0, 2).reshape(M.shape[1], -1)
# per-gene dynamic class from the profile pooled over the 4 datasets (transparent rules)
pooled = F.reshape(-1, 4, 4).mean(1)            # g x (EV, ERF5, CRF4, both)
def classify(m):
    e5, c4, bo = m[1], m[2], m[3]
    top_single = max(abs(e5), abs(c4))
    if abs(bo) > top_single + 0.25 and (abs(e5) < 0.35 and abs(c4) < 0.35 or np.sign(bo) != np.sign(e5 + c4)):
        return "combination-specific " + ("induction" if bo > 0 else "repression")
    if abs(bo) > top_single + 0.15:
        return "combination-enhanced " + ("induction" if bo > 0 else "repression")
    if abs(c4) >= abs(e5):
        return "CRF4-driven " + ("induction" if c4 > 0 else "repression")
    return "ERF5-driven " + ("induction" if e5 > 0 else "repression")
sel["module_label"] = [classify(m) for m in pooled]
codes = {l: i + 1 for i, l in enumerate(sorted(sel.module_label.unique(), key=lambda x: (not x.startswith("comb"), x)))}
sel["module"] = sel.module_label.map(codes)
sel["module_name"] = [f"M{m}: {l} (n={int((sel.module == m).sum())})" for m, l in zip(sel.module, sel.module_label)]
print(sel.module_name.value_counts().to_string())
res.to_csv(os.path.join(OUT, "concordance_all_genes.tsv"), sep="\t", index=False)
sel.to_csv(os.path.join(OUT, "concordant_genes_modules.tsv"), sep="\t", index=False)
long = []
for i, n in enumerate(names):
    for j, s in enumerate(["EV", "ERF5", "CRF4", "both"]):
        long.append(pd.DataFrame({"gene": sel.gene.values, "module_name": sel.module_name.values, "dataset": n, "state": s, "value": M[i, :, j]}))
pd.concat(long).to_csv(os.path.join(OUT, "module_profiles_long.tsv"), sep="\t", index=False)
pd.DataFrame({"null_score": np.quantile(null, np.linspace(0, 1, 2001))}).to_csv(os.path.join(OUT, "null_quantiles.tsv"), sep="\t", index=False)
