#!/usr/bin/env python3
"""Step 5: main effects + universal TF x TF interaction screen in each
experiment, fixed-effect meta-analysis across experiments, replication
statistics, and a stratified permutation null for every count reported.

Model (per gene, per experiment; log1p CP10k; HC1-robust SEs):
  y = a + sum_t b_t x_t + Z c                              (main effects)
  y = a + sum_t b_t x_t + g_ab x_a x_b + Z c               (one model per pair a,b)
x_t = GR-TF presence (>=1 junction read, the threshold that maximised
cross-experiment reproducibility in step 3); Z = log UMI, log genes detected,
log RFP, log total GR-TF reads, library (Exp2) and 9 marker-based cell-type
scores.  Interaction coefficients use Frisch-Waugh-Lovell so all 253 pairs x
all genes are fitted exactly and fast.

Consensus = meta-analysis BH q < 0.05 AND p < 0.05 with the same sign in each
experiment separately.  Empirical FDR of consensus sets from N_PERM
permutations of TF profiles within library x GR-TF-load strata.
"""
import os, sys, time, itertools, json
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from common import *

SCHEME = os.environ.get("SCHEME", "ge1")
N_PERM = int(os.environ.get("N_PERM", 20))
OUT = os.path.join(RES, f"consensus_{SCHEME}" + ("" if os.environ.get("FLEX", "1") == "1" else "_linearload")); os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(1234)

ds = {e: load(e) for e in EXPS}
genes, tfs = analysis_sets(ds)
Y = {e: response(ds[e], genes) for e in EXPS}
CT = {e: pd.read_csv(os.path.join(RES, f"{e}_celltype_scores.tsv"), sep="\t").drop(columns=["cell_type", "cell"]) for e in EXPS}
FLEX = os.environ.get("FLEX", "1") == "1"
Z = {e: np.column_stack([covariates(ds[e], tfs, flexible_load=FLEX).values, CT[e].values]) for e in EXPS}
X = {e: encode(ds[e]["T"], tfs, SCHEME) for e in EXPS}
pairs = list(itertools.combinations(range(len(tfs)), 2))
pair_names = [f"{tfs[a]}-{tfs[b]}" for a, b in pairs]


MIN_CELLS = int(os.environ.get("MIN_CELLS", 30))   # min cells in ++, +-, -+ groups


def fit_exp(Xe, Ze, Ye):
    m = ols(Xe, Ze, Ye)
    D = m["design"]; E = m["resid"]; n, p = D.shape
    W = np.column_stack([Xe[:, a] * Xe[:, b] for a, b in pairs])
    Wt = W - D @ np.linalg.lstsq(D, W, rcond=None)[0]
    ww = (Wt ** 2).sum(0)                                   # (P,)
    # support: enough cells in each informative group
    xb = Xe > 0
    sup = np.array([min((xb[:, a] & xb[:, b]).sum(), (xb[:, a] & ~xb[:, b]).sum(),
                        (~xb[:, a] & xb[:, b]).sum()) for a, b in pairs])
    ok = (ww > 1e-8) & (sup >= MIN_CELLS)
    g = (Wt.T @ E) / np.where(ok, ww, np.nan)[:, None]
    # HC3: leverage of the full model = base leverage + w~^2 / w~'w~
    h = m["lev"][:, None] + Wt ** 2 / ww[None, :]
    U = Wt ** 2 / (1 - np.minimum(h, 0.99)) ** 2            # n x P
    S2 = U.T @ (E ** 2)
    S3 = (U * Wt).T @ E
    S4 = (U * Wt ** 2).sum(0)
    var = (S2 - 2 * g * S3 + g ** 2 * S4[:, None]) / ww[:, None] ** 2
    se = np.sqrt(np.maximum(var, 1e-300))
    t = g / se
    pv = 2 * stats.t.sf(np.abs(t), n - p - 1)
    pv[~ok] = np.nan
    return dict(main=dict(beta=m["beta"], se=m["se"], p=m["p"], t=m["t"]),
                inter=dict(beta=g, se=se, p=pv, t=t), support=sup)


def meta(r1, r2):
    w1, w2 = 1 / r1["se"] ** 2, 1 / r2["se"] ** 2
    b = (r1["beta"] * w1 + r2["beta"] * w2) / (w1 + w2)
    se = np.sqrt(1 / (w1 + w2))
    z = b / se
    p = 2 * stats.norm.sf(np.abs(z))
    Q = w1 * (r1["beta"] - b) ** 2 + w2 * (r2["beta"] - b) ** 2
    return dict(beta=b, se=se, z=z, p=p, q=bh(p), p_het=stats.chi2.sf(Q, 1))


def consensus(r1, r2, mt):
    same = np.sign(r1["beta"]) == np.sign(r2["beta"])
    return (mt["q"] < 0.05) & (r1["p"] < 0.05) & (r2["p"] < 0.05) & same


def pi1(p):
    """Storey pi1 = 1 - pi0 (lambda = 0.5)."""
    p = np.asarray(p).ravel()
    return 1 - min(1, np.mean(p > 0.5) / 0.5) if p.size else np.nan


def replication(r1, r2):
    out = {}
    for name, a, b in (("1->2", r1, r2), ("2->1", r2, r1)):
        disc = bh(a["p"]) < 0.05
        same = np.sign(a["beta"]) == np.sign(b["beta"])
        out[f"n_disc_{name}"] = int(disc.sum())
        out[f"rep_rate_{name}"] = float(np.mean((b["p"][disc] < 0.05) & same[disc])) if disc.any() else np.nan
        out[f"sign_agree_{name}"] = float(np.mean(same[disc])) if disc.any() else np.nan
        # pi1 of one-sided (same direction) p-values in the other experiment
        p1s = np.where(same, b["p"] / 2, 1 - b["p"] / 2)[disc]
        out[f"pi1_{name}"] = pi1(p1s)
    return out


def summarise(F, label):
    s = {}
    for kind in ("main", "inter"):
        r1, r2 = F["Exp1"][kind], F["Exp2"][kind]
        mt = meta(r1, r2)
        cons = consensus(r1, r2, mt)
        s[kind] = dict(n_sig_Exp1=int((bh(r1["p"]) < .05).sum()), n_sig_Exp2=int((bh(r2["p"]) < .05).sum()),
                       n_meta_q05=int((mt["q"] < .05).sum()), n_consensus=int(cons.sum()),
                       r_beta=float(pd.Series(r1["t"].ravel()).corr(pd.Series(r2["t"].ravel()))),
                       **replication(r1, r2))
    return s


t0 = time.time()
F = {e: fit_exp(X[e], Z[e], Y[e]) for e in EXPS}
real = summarise(F, "real")
print("REAL", json.dumps(real, indent=1), f"[{time.time()-t0:.0f}s]", flush=True)

# ------------------------------------------------------------ save tables
def long_table(kind, names):
    r1, r2 = F["Exp1"][kind], F["Exp2"][kind]
    mt = meta(r1, r2)
    cons = consensus(r1, r2, mt)
    keep = (mt["q"] < 0.05) | (bh(r1["p"]) < 0.05) | (bh(r2["p"]) < 0.05)
    ii, jj = np.where(keep)
    df = pd.DataFrame({
        "regulator": np.array(names)[ii], "target": np.array(genes)[jj],
        "beta_Exp1": r1["beta"][ii, jj], "p_Exp1": r1["p"][ii, jj], "q_Exp1": bh(r1["p"])[ii, jj],
        "beta_Exp2": r2["beta"][ii, jj], "p_Exp2": r2["p"][ii, jj], "q_Exp2": bh(r2["p"])[ii, jj],
        "beta_meta": mt["beta"][ii, jj], "se_meta": mt["se"][ii, jj], "p_meta": mt["p"][ii, jj],
        "q_meta": mt["q"][ii, jj], "p_heterogeneity": mt["p_het"][ii, jj], "consensus": cons[ii, jj]})
    return df.sort_values(["consensus", "p_meta"], ascending=[False, True]), cons, mt


main_df, main_cons, main_meta = long_table("main", tfs)
int_df, int_cons, int_meta = long_table("inter", pair_names)
main_df.to_csv(os.path.join(OUT, "main_effects_all_significant.tsv.gz"), sep="\t", index=False)
int_df.to_csv(os.path.join(OUT, "interactions_all_significant.tsv.gz"), sep="\t", index=False)
np.savez_compressed(os.path.join(OUT, "full_stats.npz"),
                    **{f"{e}_{k}_{s}": F[e][k][s] for e in EXPS for k in ("main", "inter") for s in ("beta", "se", "p")},
                    genes=np.array(genes), tfs=np.array(tfs), pairs=np.array(pair_names))

# per-TF / per-pair summaries incl. number of supporting cells
rows = []
for j, tf in enumerate(tfs):
    r = dict(tf=tf)
    for e in EXPS:
        r[f"cells_pos_{e}"] = int(X[e][:, j].sum()) if SCHEME.startswith("ge") else np.nan
        r[f"n_q05_{e}"] = int((bh(F[e]["main"]["p"])[j] < .05).sum())
    r["r_t_Exp1_Exp2"] = float(np.corrcoef(F["Exp1"]["main"]["t"][j], F["Exp2"]["main"]["t"][j])[0, 1])
    r["n_consensus"] = int(main_cons[j].sum())
    r["n_consensus_up"] = int((main_cons[j] & (main_meta["beta"][j] > 0)).sum())
    r["n_consensus_down"] = int((main_cons[j] & (main_meta["beta"][j] < 0)).sum())
    rows.append(r)
pd.DataFrame(rows).to_csv(os.path.join(OUT, "per_TF_summary.tsv"), sep="\t", index=False)

rows = []
for k, (a, b) in enumerate(pairs):
    r = dict(pair=pair_names[k], tf_a=tfs[a], tf_b=tfs[b])
    for e in EXPS:
        xa, xb = X[e][:, a] > 0, X[e][:, b] > 0
        r.update({f"n_++_{e}": int((xa & xb).sum()), f"n_+-_{e}": int((xa & ~xb).sum()),
                  f"n_-+_{e}": int((~xa & xb).sum()), f"n_--_{e}": int((~xa & ~xb).sum()),
                  f"n_inter_q05_{e}": int((bh(F[e]["inter"]["p"])[k] < .05).sum())})
    r["n_consensus_interaction"] = int(int_cons[k].sum())
    r["n_synergy(+)"] = int((int_cons[k] & (int_meta["beta"][k] > 0)).sum())
    r["n_antagonism(-)"] = int((int_cons[k] & (int_meta["beta"][k] < 0)).sum())
    # combinatorial (additive) co-targets: consensus main-effect target of both TFs
    r["n_consensus_cotargets"] = int((main_cons[a] & main_cons[b]).sum())
    rows.append(r)
pd.DataFrame(rows).sort_values("n_consensus_interaction", ascending=False).to_csv(
    os.path.join(OUT, "per_pair_summary.tsv"), sep="\t", index=False)

# ------------------------------------------------------------ permutations
null = []
for i in range(N_PERM):
    t1 = time.time()
    Fp = {e: fit_exp(X[e][strata_permutation(ds[e], rng)], Z[e], Y[e]) for e in EXPS}
    s = summarise(Fp, "perm")
    null.append({f"{k}_{m}": v for k in s for m, v in s[k].items()})
    print(f"perm {i}: main consensus {s['main']['n_consensus']}, inter consensus {s['inter']['n_consensus']}, "
          f"main sig1/2 {s['main']['n_sig_Exp1']}/{s['main']['n_sig_Exp2']}, inter sig1/2 "
          f"{s['inter']['n_sig_Exp1']}/{s['inter']['n_sig_Exp2']} [{time.time()-t1:.0f}s]", flush=True)
    pd.DataFrame(null).to_csv(os.path.join(OUT, "permutation_null.tsv"), sep="\t", index=False)

nd = pd.DataFrame(null)
summary = dict(scheme=SCHEME, n_genes=len(genes), tfs=tfs, n_pairs=len(pairs), real=real,
               null_mean={c: float(nd[c].mean()) for c in nd.columns} if len(nd) else {},
               null_max={c: float(nd[c].max()) for c in nd.columns} if len(nd) else {})
for kind in ("main", "inter"):
    if len(nd):
        summary[f"empirical_FDR_consensus_{kind}"] = float(nd[f"{kind}_n_consensus"].mean() / max(real[kind]["n_consensus"], 1))
json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
print(json.dumps({k: v for k, v in summary.items() if k.startswith("empirical")}, indent=1))
