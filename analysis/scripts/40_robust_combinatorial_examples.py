#!/usr/bin/env python3
"""Step 40: example genes for the robust dose-based TF x TF interactions (log_loadmod,
strict tier, surviving non-linear dose terms, pair profile replicated genome-wide).
Exports 4-state values (neither / TFa / TFb / both; presence = >=1 read, for display)
per experiment as 10 random pools of cells, adjusted for depth, load (incl. per-TF
load modulation), cell type and the other TFs."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
os.environ.setdefault("N_PERM", "0")
src = open(os.path.join(os.path.dirname(__file__), "05_consensus.py")).read()
exec(compile(src[:src.index("t0 = time.time()")], "05", "exec"))
OUT = os.path.join(RES, "combinatorial_push"); rng = np.random.default_rng(40)
calls = pd.read_csv(os.path.join(OUT, "calls_log_loadmod.tsv.gz"), sep="\t")
pr = pd.read_csv(os.path.join(OUT, "pair_replication.tsv"), sep="\t"); pr = pr[pr.setting == "log_loadmod"].set_index("pair")
nl = pd.read_csv(os.path.join(OUT, "nonlinear_dose_check.tsv"), sep="\t").set_index("pair")
c = calls[calls.tier == "strict"].copy()
c["pair_q"] = c.pair.map(pr.q)
c["score"] = -np.log10(np.maximum(c.p_Exp1, c.p_Exp2))
c = c[c.pair_q < 0.05].sort_values("score", ascending=False)
sel = c.groupby("pair").head(2).head(12)
print(sel[["pair", "target", "beta_Exp1", "beta_Exp2", "p_Exp1", "p_Exp2", "pair_q"]].round(4).to_string(index=False))
sel.to_csv(os.path.join(OUT, "examples_selected.tsv"), sep="\t", index=False)
Xd = {e: encode(ds[e]["T"], tfs, "log") for e in EXPS}
L = {e: covariates(ds[e], tfs)["log_tf_total"].values for e in EXPS}; L = {e: (v - v.mean()) / v.std() for e, v in L.items()}
gi = {g: i for i, g in enumerate(genes)}; tfi = {t: i for i, t in enumerate(tfs)}
rows = []
for _, r in sel.iterrows():
    a, b = (tfi[x] for x in r.pair.split("-")); g = gi[r.target]
    for e in EXPS:
        xa, xb = (Xd[e][:, a] > 0).astype(float), (Xd[e][:, b] > 0).astype(float)
        S = np.column_stack([xa * (1 - xb), (1 - xa) * xb, xa * xb])
        others = np.delete(Xd[e], [a, b], 1)
        D = np.column_stack([np.ones(len(xa)), S, others, Z[e], Xd[e] * L[e][:, None]])
        y = Y[e][:, g]; B = np.linalg.lstsq(D, y, rcond=None)[0]
        adj = y - D[:, 4:] @ B[4:]
        state = np.select([xa * xb == 1, xa == 1, xb == 1], ["both", "a", "b"], "neither")
        ref = adj[state == "neither"].mean()
        for st in ("neither", "a", "b", "both"):
            for k, p in enumerate(np.array_split(rng.permutation(np.where(state == st)[0]), 10)):
                rows.append(dict(pair=r.pair, TFa=r.pair.split("-")[0], TFb=r.pair.split("-")[1], target=r.target,
                                 exp=f"Single-cell {e}", state=st, rel=adj[p].mean() - ref, n=int((state == st).sum())))
pd.DataFrame(rows).to_csv(os.path.join(OUT, "examples_values.tsv"), sep="\t", index=False)
