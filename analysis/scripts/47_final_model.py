#!/usr/bin/env python3
"""Step 47: FINAL model for the revision (raw counts, QC cells, both experiments).
TF encoding  : dose = log1p(GR-TF junction reads)
Covariates   : log UMI, log genes, log RFP, log total GR-TF reads (+ squared/cubic
               terms and their product), library, 9 marker cell-type scores, and
               each TF's dose x standardised GR-TF load ("load-modulated" baseline)
Main effects : y ~ 23 TF doses + covariates          (per experiment, HC3 SEs)
Interactions : + one TFa x TFb dose product at a time (253 pairs, FWL, HC3)
Combination  : fixed-effect meta-analysis; tiers
               strict : meta BH q<0.05 and p<0.05, same sign in each experiment
               medium : meta BH q<0.10 and p<0.10, same sign in each experiment
Null         : N_PERM permutations of TF profiles within library x load-decile strata
"""
import os, sys, json, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
os.environ["N_PERM"] = "0"; os.environ["SCHEME"] = "log"
src = open(os.path.join(os.path.dirname(__file__), "05_consensus.py")).read()
exec(compile(src[:src.index("t0 = time.time()")], "05", "exec"))
NP = int(os.environ.get("NP", 50))
OUT = os.path.join(RES, "final"); os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(2026)
L = {e: covariates(ds[e], tfs)["log_tf_total"].values for e in EXPS}; L = {e: (v - v.mean()) / v.std() for e, v in L.items()}
TIERS = {"strict": (0.05, 0.05), "medium": (0.10, 0.10)}


def fit_all(perm=False):
    F = {}
    for e in EXPS:
        Xe = X[e][strata_permutation(ds[e], rng)] if perm else X[e]
        F[e] = fit_exp(Xe, np.column_stack([Z[e], Xe * L[e][:, None]]), Y[e])
    return F


def tiers(r1, r2, mt):
    same = np.sign(r1["beta"]) == np.sign(r2["beta"])
    return {t: (mt["q"] < q) & (r1["p"] < p) & (r2["p"] < p) & same for t, (q, p) in TIERS.items()}


def counts(F):
    out = {}
    for k in ("main", "inter"):
        r1, r2 = F["Exp1"][k], F["Exp2"][k]; mt = meta(r1, r2)
        for t, c in tiers(r1, r2, mt).items():
            out[f"{k}_{t}"] = int(c.sum())
    return out


t0 = time.time()
F = fit_all()
real = counts(F); print("REAL", real, f"[{time.time()-t0:.0f}s]", flush=True)
null = []
for i in range(NP):
    null.append(counts(fit_all(True)))
    if i % 10 == 9:
        print("perm", i + 1, pd.DataFrame(null).mean().round(2).to_dict(), flush=True)
nd = pd.DataFrame(null); nd.to_csv(os.path.join(OUT, "permutation_null.tsv"), sep="\t", index=False)
summary = dict(real=real, null_mean=nd.mean().to_dict(), null_max=nd.max().to_dict(),
               FDR={k: float(nd[k].mean() / max(real[k], 1)) for k in real}, n_perm=NP,
               cells={e: int(len(Y[e])) for e in EXPS}, genes=len(genes), tfs=tfs)
json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=1)

# ---- tables
for k, names in (("main", tfs), ("inter", pair_names)):
    r1, r2 = F["Exp1"][k], F["Exp2"][k]; mt = meta(r1, r2); T = tiers(r1, r2, mt)
    keep = T["medium"] | (mt["q"] < 0.1)
    ii, jj = np.where(keep)
    df = pd.DataFrame({"regulator" if k == "main" else "pair": np.array(names)[ii], "target": np.array(genes)[jj],
                       "beta_Exp1": r1["beta"][ii, jj], "p_Exp1": r1["p"][ii, jj], "beta_Exp2": r2["beta"][ii, jj], "p_Exp2": r2["p"][ii, jj],
                       "beta_meta": mt["beta"][ii, jj], "q_meta": mt["q"][ii, jj], "p_heterogeneity": mt["p_het"][ii, jj],
                       "tier": np.where(T["strict"][ii, jj], "strict", np.where(T["medium"][ii, jj], "medium", "meta q<0.1 only"))})
    df.sort_values(["tier", "q_meta"]).to_csv(os.path.join(OUT, f"{k}_effects.tsv.gz"), sep="\t", index=False)
np.savez_compressed(os.path.join(OUT, "full_stats.npz"), **{f"{e}_{k}_{s}": F[e][k][s] for e in EXPS for k in ("main", "inter") for s in ("beta", "se", "p", "t")},
                    genes=np.array(genes), tfs=np.array(tfs), pairs=np.array(pair_names))
# per TF / per pair summaries
r1, r2 = F["Exp1"]["main"], F["Exp2"]["main"]; mt = meta(r1, r2); T = tiers(r1, r2, mt)
rows = []
for j, tf in enumerate(tfs):
    rows.append(dict(tf=tf, cells_pos_Exp1=int((X["Exp1"][:, j] > 0).sum()), cells_pos_Exp2=int((X["Exp2"][:, j] > 0).sum()),
                     frac_pos_Exp1=float((X["Exp1"][:, j] > 0).mean()), frac_pos_Exp2=float((X["Exp2"][:, j] > 0).mean()),
                     r_Exp1_Exp2=float(np.corrcoef(r1["t"][j], r2["t"][j])[0, 1]),
                     consensus_strict=int(T["strict"][j].sum()), up=int((T["strict"][j] & (mt["beta"][j] > 0)).sum()),
                     down=int((T["strict"][j] & (mt["beta"][j] < 0)).sum())))
pd.DataFrame(rows).sort_values("consensus_strict", ascending=False).to_csv(os.path.join(OUT, "per_TF.tsv"), sep="\t", index=False)
r1, r2 = F["Exp1"]["inter"], F["Exp2"]["inter"]; mt = meta(r1, r2); T = tiers(r1, r2, mt)
rows = []
for k, (a, b) in enumerate(pairs):
    xa1, xb1, xa2, xb2 = X["Exp1"][:, a] > 0, X["Exp1"][:, b] > 0, X["Exp2"][:, a] > 0, X["Exp2"][:, b] > 0
    rows.append(dict(pair=pair_names[k], n_both_Exp1=int((xa1 & xb1).sum()), n_a_only_Exp1=int((xa1 & ~xb1).sum()), n_b_only_Exp1=int((~xa1 & xb1).sum()), n_neither_Exp1=int((~xa1 & ~xb1).sum()),
                     n_both_Exp2=int((xa2 & xb2).sum()), n_a_only_Exp2=int((xa2 & ~xb2).sum()), n_b_only_Exp2=int((~xa2 & xb2).sum()), n_neither_Exp2=int((~xa2 & ~xb2).sum()),
                     interactions_strict=int(T["strict"][k].sum()), interactions_medium=int(T["medium"][k].sum())))
pd.DataFrame(rows).sort_values("interactions_medium", ascending=False).to_csv(os.path.join(OUT, "per_pair.tsv"), sep="\t", index=False)
print(json.dumps(summary["FDR"], indent=1))
