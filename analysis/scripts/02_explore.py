#!/usr/bin/env python3
"""Step 2: describe GR-TF detection in both experiments (counts, load, RFP,
co-occurrence, endogenous-locus check). Prints a report and writes tables."""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from common import *

rows = []
for exp in EXPS:
    d = load(exp)
    G, T, rfp = d["G"], d["T"], d["rfp"]
    lib = G.sum(0).A1
    ngenes = (G > 0).sum(0).A1
    tot_tf = T.sum(1).values
    print(f"\n===== {exp}: {T.shape[0]} cells")
    print("cor(log total TF reads, log gene UMI) = %.2f" % np.corrcoef(np.log1p(tot_tf), np.log(lib))[0, 1])
    print("cor(log total TF reads, log RFP)      = %.2f" % np.corrcoef(np.log1p(tot_tf), np.log1p(rfp))[0, 1])
    print("cor(log RFP, log gene UMI)            = %.2f" % np.corrcoef(np.log1p(rfp), np.log(lib))[0, 1])
    print("RFP quantiles", np.quantile(rfp, [0, .05, .5, .95, 1]).round(0))
    print("total TF reads quantiles", np.quantile(tot_tf, [0, .05, .25, .5, .75, .95, 1]).round(0))
    for k in (1, 2, 3, 5):
        n = (T >= k).sum(1)
        print(f"#TFs with >= {k} reads per cell: median {np.median(n):.0f}, 0-TF cells {np.mean(n==0)*100:.1f}%")
    # per TF stats
    gidx = {g: i for i, g in enumerate(d["genes"])}
    for tf in T.columns:
        v = T[tf].values
        r = dict(exp=exp, tf=tf, total=int(v.sum()),
                 frac_ge1=np.mean(v >= 1), frac_ge2=np.mean(v >= 2), frac_ge3=np.mean(v >= 3),
                 frac_ge5=np.mean(v >= 5), median_pos=np.median(v[v > 0]) if (v > 0).any() else 0,
                 share_of_TF_reads=v.sum() / tot_tf.sum())
        # share of each cell's TF reads (composition)
        with np.errstate(invalid="ignore", divide="ignore"):
            r["mean_share_in_pos_cells"] = np.nanmean(np.where(v > 0, v / tot_tf, np.nan))
        agi = TF_AGI.get(tf)
        if agi in gidx:
            e = G[gidx[agi], :].toarray().ravel()
            r["endog_locus_frac_cells"] = np.mean(e > 0)
            r["endog_vs_tf_spearman"] = stats.spearmanr(e, v).correlation
        rows.append(r)
    # co-occurrence vs independence given load: observed pair co-detection
    # vs expectation from a model where each TF is detected independently
    # with probability depending on the cell's total TF reads (logistic fit).
    P = (T.values >= 1).astype(float)
    x = np.log1p(tot_tf)
    exp_pairs = []
    pred = np.zeros_like(P)
    for j in range(P.shape[1]):
        if P[:, j].sum() < 10:
            pred[:, j] = P[:, j].mean(); continue
        import statsmodels.api as sm
        m = sm.GLM(P[:, j], sm.add_constant(x), family=sm.families.Binomial()).fit()
        pred[:, j] = m.predict(sm.add_constant(x))
    obs = P.T @ P; ex = pred.T @ pred
    iu = np.triu_indices(P.shape[1], 1)
    ratio = (obs[iu] + 1) / (ex[iu] + 1)
    print("pairwise co-detection obs/expected (given load): median %.3f, 5-95%% [%.3f, %.3f]" %
          (np.median(ratio), *np.quantile(ratio, [.05, .95])))
    print("raw phi correlations between TF presences: median %.3f" %
          np.median(np.corrcoef(P.T)[iu]))
    # residual correlation after load
    R = P - pred
    rc = np.corrcoef(R.T)[iu]
    print("residual (load-adjusted) presence correlations: median %.3f, 5-95%% [%.3f, %.3f]" %
          (np.nanmedian(rc), *np.nanquantile(rc, [.05, .95])))

tab = pd.DataFrame(rows)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(tab.round(3).to_string())
tab.to_csv(os.path.join(RES, "tf_detection_summary.tsv"), sep="\t", index=False)
