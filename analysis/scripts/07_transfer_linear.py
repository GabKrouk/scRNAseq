#!/usr/bin/env python3
"""Step 7: cross-experiment prediction (train on one experiment, predict the
other) for linear models, with the ridge penalty calibrated on transfer.

Models (features = GR-TF presence >=1 read; covariates Z partialled out of
features and expression within each experiment, Frisch-Waugh-Lovell):
  additive      : 23 TF main effects (ridge, penalty scanned)
  interactions  : best additive model + 253 TF x TF products orthogonalised to
                  the main effects (ridge on the additive residuals, own penalty)
For each direction (Exp1->Exp2, Exp2->Exp1) and each penalty lambda we fit
on all cells of the training experiment and report, per gene, the correlation
between predicted and observed expression in the test experiment.  A
permutation of the test experiment's TF profiles gives the null.
"""
import os, sys, itertools, json, time
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from common import *

rng = np.random.default_rng(99)
ds = {e: load(e) for e in EXPS}
genes, tfs = analysis_sets(ds)
pairs = list(itertools.combinations(range(len(tfs)), 2))
CT = {e: pd.read_csv(os.path.join(RES, f"{e}_celltype_scores.tsv"), sep="\t").drop(columns=["cell_type", "cell"]) for e in EXPS}
Z = {e: np.column_stack([np.ones(ds[e]["T"].shape[0]), covariates(ds[e], tfs, flexible_load=True).values, CT[e].values]) for e in EXPS}


def resid(A, D):
    return A - D @ np.linalg.lstsq(D, A, rcond=None)[0]


def features(e, perm=None):
    """[X main effects | X_a*X_b products orthogonalised to the main effects],
    both partialled on Z: the product block then carries ONLY non-additive
    (interaction) information, so 'interactions' can only add to 'additive'
    if TF x TF interactions reproduce across experiments."""
    X = encode(ds[e]["T"], tfs, os.environ.get("SCHEME", "ge1"))
    if perm is not None:
        X = X[perm]
    W = np.column_stack([X[:, a] * X[:, b] for a, b in pairs])
    Xr = resid(X, Z[e])
    Wr = resid(W, np.column_stack([Z[e], X]))
    return np.column_stack([Xr, Wr])


Yt = {e: resid(response(ds[e], genes), Z[e]) for e in EXPS}
F = {e: features(e) for e in EXPS}
LAMBDAS = [1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1e3, 1e4, 1e6]


def corr_cols(P, Y):
    P = P - P.mean(0); Yc = Y - Y.mean(0)
    den = np.sqrt((P ** 2).sum(0) * (Yc ** 2).sum(0))
    return np.where(den > 0, (P * Yc).sum(0) / np.where(den > 0, den, 1), 0)


def ridge_pred(A, B, Ytr, lam_list):
    sd = A.std(0); ok = sd > 1e-6
    A = A[:, ok] / sd[ok]; B = B[:, ok] / sd[ok]
    U, S, Vt = np.linalg.svd(A, full_matrices=False)
    M = U.T @ Ytr
    n = A.shape[0]
    for lam in lam_list:
        d = S / (S ** 2 + lam * n)
        yield lam, B @ (Vt.T @ (d[:, None] * M))


def transfer(tr, te, Fte, cols):
    """additive: ridge on main-effect block.  interactions: best additive
    prediction + ridge on the (orthogonalised) interaction block fitted to the
    additive model's training residuals, scanned over its own penalty."""
    k = len(tfs)
    add = {lam: P for lam, P in ridge_pred(F[tr][:, :k], Fte[:, :k], Yt[tr], LAMBDAS)}
    if cols == "additive":
        return {lam: corr_cols(P, Yt[te]) for lam, P in add.items()}
    lam_a = best_add_lambda[tr]
    # training residuals of the additive model at its best penalty
    Ptr_add = next(ridge_pred(F[tr][:, :k], F[tr][:, :k], Yt[tr], [lam_a]))[1]
    out = {}
    for lam, P in ridge_pred(F[tr][:, k:], Fte[:, k:], Yt[tr] - Ptr_add, LAMBDAS):
        out[lam] = corr_cols(add[lam_a] + P, Yt[te])
    return out


best_add_lambda = {}
rows, best = [], {}
main_cols = list(range(len(tfs)))
full_cols = list(range(len(tfs) + len(pairs)))
for tr, te in (("Exp1", "Exp2"), ("Exp2", "Exp1")):
    n_te = Yt[te].shape[0]
    for model, cols in (("additive", "additive"), ("interactions", "interactions")):
        t0 = time.time()
        R = transfer(tr, te, F[te], cols)
        # null: permute the test experiment's TF profiles (strata)
        Rn = transfer(tr, te, features(te, strata_permutation(ds[te], rng)), cols)
        for lam in LAMBDAS:
            r = R[lam]; rn = Rn[lam]
            z = np.arctanh(np.clip(r, -0.999, 0.999)) * np.sqrt(n_te - 3)
            p = stats.norm.sf(z)
            zn = np.arctanh(np.clip(rn, -0.999, 0.999)) * np.sqrt(n_te - 3)
            rows.append(dict(direction=f"{tr}->{te}", model=model, lam=lam, mean_r=r.mean(), median_r=np.median(r),
                             n_sig_bh=int((bh(p) < 0.05).sum()), null_mean_r=rn.mean(),
                             null_n_sig_bh=int((bh(stats.norm.sf(zn)) < 0.05).sum())))
        k = max(LAMBDAS, key=lambda l: R[l].mean())
        best[(tr, model)] = (k, R[k], Rn[k])
        if model == "additive":
            best_add_lambda[tr] = k
        print(f"{tr}->{te} {model:12s} best lambda {k}: mean r {R[k].mean():.4f} (null {Rn[k].mean():.4f}) "
              f"[{time.time()-t0:.0f}s]", flush=True)

tab = pd.DataFrame(rows)
tab.to_csv(os.path.join(RES, "transfer_linear_lambda_grid.tsv"), sep="\t", index=False)
print(tab.round(4).to_string())
per_gene = pd.DataFrame({"gene": genes})
for (tr, model), (lam, r, rn) in best.items():
    per_gene[f"r_{model}_{tr}->{'Exp2' if tr == 'Exp1' else 'Exp1'}"] = r
    per_gene[f"null_r_{model}_{tr}->{'Exp2' if tr == 'Exp1' else 'Exp1'}"] = rn
per_gene.to_csv(os.path.join(RES, "transfer_linear_per_gene.tsv.gz"), sep="\t", index=False)
json.dump({f"{k[0]}_{k[1]}": v[0] for k, v in best.items()}, open(os.path.join(RES, "transfer_linear_best_lambda.json"), "w"))
