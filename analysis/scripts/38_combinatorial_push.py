#!/usr/bin/env python3
"""Step 38: push the search for reproducible TF x TF combinatorial effects on raw
counts, with every relaxation calibrated by the same permutation null.

Grid: TF encoding (>=1 read, >=2 reads, log1p dose) x baseline model (standard, or
'load-modulated' where every TF main effect may scale with total GR-TF load).
For each setting and for N_PERM label permutations:
  * gene-level calls at three stringency tiers
      strict : meta BH q<0.05, p<0.05 & same sign in each experiment
      medium : meta BH q<0.10, p<0.10 & same sign in each
      relaxed: meta BH q<0.20, p<0.20 & same sign in each
    -> observed count, null mean, empirical FDR
  * pair-level replication: Spearman correlation, across genes, of the pair's
    interaction t-statistics in Exp1 vs Exp2 (does the whole interaction profile
    of a pair reproduce, even if single genes are weak?) -> per-pair empirical p
    from the pooled permutation null, BH across pairs.
"""
import os, sys, json, itertools
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
os.environ.setdefault("N_PERM", "0")
src = open(os.path.join(os.path.dirname(__file__), "05_consensus.py")).read()
exec(compile(src[:src.index("t0 = time.time()")], "05", "exec"))
N_PERM = int(os.environ.get("NP", 10))
OUT = os.path.join(RES, "combinatorial_push"); os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(38)
TIERS = {"strict": (0.05, 0.05), "medium": (0.10, 0.10), "relaxed": (0.20, 0.20)}
load = {e: covariates(ds[e], tfs)["log_tf_total"].values for e in EXPS}
load = {e: (v - v.mean()) / v.std() for e, v in load.items()}


def calls(F):
    r1, r2 = F["Exp1"]["inter"], F["Exp2"]["inter"]
    mt = meta(r1, r2); same = np.sign(r1["beta"]) == np.sign(r2["beta"])
    return {t: (mt["q"] < q) & (r1["p"] < p) & (r2["p"] < p) & same for t, (q, p) in TIERS.items()}, mt


def pair_corr(F):
    t1, t2 = F["Exp1"]["inter"]["t"], F["Exp2"]["inter"]["t"]
    out = np.full(len(pairs), np.nan)
    for k in range(len(pairs)):
        ok = np.isfinite(t1[k]) & np.isfinite(t2[k])
        if ok.sum() > 100:
            out[k] = stats.spearmanr(t1[k][ok], t2[k][ok]).correlation
    return out


def fit_all(Xd, loadmod, perm=False):
    F = {}
    for e in EXPS:
        Xe = Xd[e]
        if perm:
            Xe = Xe[strata_permutation(ds[e], rng)]
        Ze = Z[e]
        if loadmod:
            Ze = np.column_stack([Ze, Xe * load[e][:, None]])
        F[e] = fit_exp(Xe, Ze, Y[e])
    return F


rows, pair_rows = [], []
for scheme in ("ge1", "ge2", "log"):
    Xd = {e: encode(ds[e]["T"], tfs, scheme) for e in EXPS}
    for loadmod in (False, True):
        tag = f"{scheme}{'_loadmod' if loadmod else ''}"
        F = fit_all(Xd, loadmod)
        C, mt = calls(F); pc = pair_corr(F)
        null_c = {t: [] for t in TIERS}; null_pc = []
        for i in range(N_PERM):
            Fp = fit_all(Xd, loadmod, perm=True)
            Cp, _ = calls(Fp)
            for t in TIERS:
                null_c[t].append(int(Cp[t].sum()))
            null_pc.append(pair_corr(Fp))
        null_pc = np.concatenate(null_pc); null_pc = null_pc[np.isfinite(null_pc)]
        pemp = np.array([(1 + np.sum(null_pc >= r)) / (1 + null_pc.size) if np.isfinite(r) else np.nan for r in pc])
        qp = bh(pemp)
        r = dict(setting=tag, n_pairs_tested=int(np.isfinite(pc).sum()),
                 pairs_q10=int((qp < 0.10).sum()), pairs_q25=int((qp < 0.25).sum()),
                 null_pair_r_95=float(np.quantile(null_pc, .95)), median_pair_r=float(np.nanmedian(pc)))
        for t in TIERS:
            n = int(C[t].sum()); nm = float(np.mean(null_c[t])) if N_PERM else np.nan
            r.update({f"{t}_n": n, f"{t}_null": nm, f"{t}_FDR": nm / max(n, 1)})
        rows.append(r)
        print(json.dumps(r), flush=True)
        for k, pr in enumerate(pairs):
            pair_rows.append(dict(setting=tag, pair=pair_names[k], r_exp1_exp2=pc[k], p_emp=pemp[k], q=qp[k],
                                  **{f"n_{t}": int(C[t][k].sum()) for t in TIERS}))
        # save gene-level calls for this setting
        ii, jj = np.where(C["relaxed"])
        pd.DataFrame({"pair": np.array(pair_names)[ii], "target": np.array(genes)[jj],
                      "beta_Exp1": F["Exp1"]["inter"]["beta"][ii, jj], "p_Exp1": F["Exp1"]["inter"]["p"][ii, jj],
                      "beta_Exp2": F["Exp2"]["inter"]["beta"][ii, jj], "p_Exp2": F["Exp2"]["inter"]["p"][ii, jj],
                      "q_meta": mt["q"][ii, jj], "tier": np.where(C["strict"][ii, jj], "strict", np.where(C["medium"][ii, jj], "medium", "relaxed"))}
                     ).to_csv(os.path.join(OUT, f"calls_{tag}.tsv.gz"), sep="\t", index=False)
pd.DataFrame(rows).to_csv(os.path.join(OUT, "grid_summary.tsv"), sep="\t", index=False)
pd.DataFrame(pair_rows).to_csv(os.path.join(OUT, "pair_replication.tsv"), sep="\t", index=False)
