#!/usr/bin/env python3
"""Step 45: values for the Sup Info box plots of the 461 reproducible non-additive
TF x TF interactions (log_loadmod, strict + medium). For each interaction and each
experiment: cells in the 4 presence states (neither / TFa only / TFb only / both,
>=1 junction read), expression adjusted for depth, load (incl. per-TF load
modulation), cell type and the other 21 TFs, summarised as 10 random pools per state
and expressed relative to the mean of cells with neither TF."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
os.environ.setdefault("N_PERM", "0")
src = open(os.path.join(os.path.dirname(__file__), "05_consensus.py")).read()
exec(compile(src[:src.index("t0 = time.time()")], "05", "exec"))
OUT = os.path.join(RES, "combinatorial_push"); rng = np.random.default_rng(45)
calls = pd.read_csv(os.path.join(OUT, "calls_log_loadmod.tsv.gz"), sep="\t")
sel = calls[calls.tier.isin(["strict", "medium"])].copy()
na = pd.read_csv(os.path.join(OUT, "nonadditivity_medium.tsv"), sep="\t")
sel = sel.merge(na.groupby(["pair", "target"]).nonadd.mean().rename("nonadd").reset_index(), on=["pair", "target"])
sel["type"] = np.where(sel.nonadd > 0, "synergy", "antagonism")
sel.to_csv(os.path.join(OUT, "supinfo_461_list.tsv"), sep="\t", index=False)
Xd = {e: encode(ds[e]["T"], tfs, "log") for e in EXPS}
L = {e: covariates(ds[e], tfs)["log_tf_total"].values for e in EXPS}; L = {e: (v - v.mean()) / v.std() for e, v in L.items()}
gi = {g: i for i, g in enumerate(genes)}; tfi = {t: i for i, t in enumerate(tfs)}
rows = []
for pair, sub in sel.groupby("pair"):
    a, b = (tfi[x] for x in pair.split("-")); gidx = [gi[g] for g in sub.target]
    for e in EXPS:
        xa, xb = (Xd[e][:, a] > 0).astype(float), (Xd[e][:, b] > 0).astype(float)
        S = np.column_stack([xa * (1 - xb), (1 - xa) * xb, xa * xb])
        D = np.column_stack([np.ones(len(xa)), S, np.delete(Xd[e], [a, b], 1), Z[e], Xd[e] * L[e][:, None]])
        Yg = Y[e][:, gidx]; B = np.linalg.lstsq(D, Yg, rcond=None)[0]
        adj = Yg - D[:, 4:] @ B[4:]
        state = np.select([xa * xb == 1, xa == 1, xb == 1], ["both", "a", "b"], "neither")
        ref = adj[state == "neither"].mean(0)
        for st in ("neither", "a", "b", "both"):
            for k, p in enumerate(np.array_split(rng.permutation(np.where(state == st)[0]), 10)):
                v = adj[p].mean(0) - ref
                for j, g in enumerate(sub.target):
                    rows.append((pair, g, e, st, k, v[j]))
df = pd.DataFrame(rows, columns=["pair", "target", "exp", "state", "pool", "rel"])
df.to_csv(os.path.join(OUT, "supinfo_461_values.tsv.gz"), sep="\t", index=False)
print(len(sel), "interactions;", df.shape)
