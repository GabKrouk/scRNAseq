#!/usr/bin/env python3
"""Step 39: are the dose-based interactions (log_loadmod setting) really pair-specific,
or non-linear single-TF dose responses? For each pair with calls, refit the
interaction adding squared and cubic dose terms of BOTH TFs, then re-apply the same
tier thresholds. Also reports the correlation between the two TFs' doses."""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
os.environ.setdefault("N_PERM", "0")
src = open(os.path.join(os.path.dirname(__file__), "05_consensus.py")).read()
exec(compile(src[:src.index("t0 = time.time()")], "05", "exec"))
OUT = os.path.join(RES, "combinatorial_push")
calls = pd.read_csv(os.path.join(OUT, "calls_log_loadmod.tsv.gz"), sep="\t")
pr = pd.read_csv(os.path.join(OUT, "pair_replication.tsv"), sep="\t"); pr = pr[pr.setting == "log_loadmod"]
Xd = {e: encode(ds[e]["T"], tfs, "log") for e in EXPS}
L = {e: covariates(ds[e], tfs)["log_tf_total"].values for e in EXPS}; L = {e: (v - v.mean()) / v.std() for e, v in L.items()}
Zb = {e: np.column_stack([Z[e], Xd[e] * L[e][:, None]]) for e in EXPS}
gi = {g: i for i, g in enumerate(genes)}; tfi = {t: i for i, t in enumerate(tfs)}
rows = []
for pair in sorted(set(calls[calls.tier.isin(["strict", "medium"])].pair) | set(pr[pr.q < 0.1].pair)):
    a, b = (tfi[x] for x in pair.split("-"))
    sub = calls[(calls.pair == pair)]
    gidx = [gi[g] for g in sub.target] if len(sub) else []
    res = {}
    for e in EXPS:
        xa, xb = Xd[e][:, a], Xd[e][:, b]
        extra = np.column_stack([xa ** 2, xb ** 2, xa ** 3, xb ** 3])
        D = np.column_stack([Xd[e], xa * xb])
        f = ols(D, np.column_stack([Zb[e], extra]), Y[e][:, gidx] if gidx else Y[e][:, :1])
        res[e] = (f["beta"][-1], f["se"][-1], f["p"][-1])
    r = dict(pair=pair, dose_corr_Exp1=stats.spearmanr(Xd["Exp1"][:, a], Xd["Exp1"][:, b]).correlation,
             dose_corr_Exp2=stats.spearmanr(Xd["Exp2"][:, a], Xd["Exp2"][:, b]).correlation,
             pair_q=float(pr.set_index("pair").q.get(pair, np.nan)))
    if gidx:
        b1, s1, p1 = res["Exp1"]; b2, s2, p2 = res["Exp2"]
        same = np.sign(b1) == np.sign(b2)
        for t, thr in (("strict", .05), ("medium", .10)):
            orig = sub.tier.isin(["strict"] if t == "strict" else ["strict", "medium"]).values
            keep = (p1 < thr) & (p2 < thr) & same
            r[f"{t}_orig"] = int(orig.sum()); r[f"{t}_kept_after_nonlinear_dose"] = int((orig & keep).sum())
        r["corr_beta_Exp1"] = float(np.corrcoef(sub.beta_Exp1, b1)[0, 1]) if len(sub) > 2 else np.nan
    rows.append(r)
df = pd.DataFrame(rows).sort_values("pair_q")
df.to_csv(os.path.join(OUT, "nonlinear_dose_check.tsv"), sep="\t", index=False)
pd.set_option("display.width", 220); print(df.round(3).to_string(index=False))
print("TOTAL strict kept %d/%d, medium kept %d/%d" % (df.strict_kept_after_nonlinear_dose.sum(), df.strict_orig.sum(),
      df.medium_kept_after_nonlinear_dose.sum(), df.medium_orig.sum()))
