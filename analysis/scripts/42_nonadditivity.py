#!/usr/bin/env python3
"""Step 42: show the robust interactions as NON-ADDITIVE effects (ANOVA-type
interaction = combined effect differs from the sum of the two single-TF effects).
For each strict interaction (log_loadmod setting) and each experiment, cells are
grouped by dose: none of both TFs; TFa high (>= median dose of TFa+ cells) without
TFb; TFb high without TFa; both high. Group means of adjusted expression (depth,
load, per-TF load modulation, cell type, other TFs) relative to 'none' give
dA, dB, dAB; additive expectation = dA + dB; non-additivity = dAB - (dA + dB).
Also exports 3x3 dose landscapes (0 / low / high for each TF) for example genes."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
os.environ.setdefault("N_PERM", "0")
src = open(os.path.join(os.path.dirname(__file__), "05_consensus.py")).read()
exec(compile(src[:src.index("t0 = time.time()")], "05", "exec"))
OUT = os.path.join(RES, "combinatorial_push")
calls = pd.read_csv(os.path.join(OUT, "calls_log_loadmod.tsv.gz"), sep="\t")
SET = os.environ.get("SET", "strict")          # strict (126) or medium (= strict + medium, 461)
strict = calls[calls.tier.isin(["strict"] if SET == "strict" else ["strict", "medium"])].reset_index(drop=True)
ex = pd.read_csv(os.path.join(OUT, "examples_selected.tsv"), sep="\t")
Xd = {e: encode(ds[e]["T"], tfs, "log") for e in EXPS}
L = {e: covariates(ds[e], tfs)["log_tf_total"].values for e in EXPS}; L = {e: (v - v.mean()) / v.std() for e, v in L.items()}
gi = {g: i for i, g in enumerate(genes)}; tfi = {t: i for i, t in enumerate(tfs)}
rows, land = [], []
for pair, sub in strict.groupby("pair"):
    a, b = (tfi[x] for x in pair.split("-"))
    gidx = [gi[g] for g in sub.target]
    for e in EXPS:
        da, db = Xd[e][:, a], Xd[e][:, b]
        others = np.delete(Xd[e], [a, b], 1)
        D = np.column_stack([np.ones(len(da)), others, Z[e], Xd[e] * L[e][:, None]])
        Yg = Y[e][:, gidx]
        adj = Yg - D[:, 1:] @ np.linalg.lstsq(D, Yg, rcond=None)[0][1:]   # remove others/covariates, keep TFa/TFb signal
        ha = da >= np.median(da[da > 0]); hb = db >= np.median(db[db > 0])
        grp = {"none": (da == 0) & (db == 0), "A": ha & (db == 0), "B": hb & (da == 0), "AB": ha & hb}
        m = {k: adj[v].mean(0) for k, v in grp.items()}
        n = {k: int(v.sum()) for k, v in grp.items()}
        for j, g in enumerate(sub.target):
            dA, dB, dAB = m["A"][j] - m["none"][j], m["B"][j] - m["none"][j], m["AB"][j] - m["none"][j]
            rows.append(dict(pair=pair, target=g, exp=e, dA=dA, dB=dB, dAB=dAB, additive=dA + dB, nonadd=dAB - dA - dB,
                             n_none=n["none"], n_A=n["A"], n_B=n["B"], n_AB=n["AB"]))
        # 3x3 dose landscapes for example genes
        for g in ex[ex.pair == pair].target:
            j = list(sub.target).index(g) if g in list(sub.target) else None
            if j is None: continue
            def lvl(d):
                pos = d[d > 0]; t = np.median(pos)
                return np.where(d == 0, "0", np.where(d < t, "low", "high"))
            la, lb = lvl(da), lvl(db)
            base = adj[(la == "0") & (lb == "0"), j].mean()
            for x in ("0", "low", "high"):
                for y in ("0", "low", "high"):
                    msk = (la == x) & (lb == y)
                    land.append(dict(pair=pair, target=g, exp=e, doseA=x, doseB=y, mean=adj[msk, j].mean() - base, n=int(msk.sum())))
df = pd.DataFrame(rows); df.to_csv(os.path.join(OUT, f"nonadditivity_{SET}.tsv"), sep="\t", index=False)
pd.DataFrame(land).to_csv(os.path.join(OUT, "dose_landscapes_examples.tsv" if SET == "strict" else f"dose_landscapes_examples_{SET}.tsv"), sep="\t", index=False)
w = df.pivot_table(index=["pair", "target"], columns="exp", values="nonadd")
print(SET, "interactions:", len(w), "| non-additivity same sign in Exp1 & Exp2: %.0f%%" % (100 * np.mean(np.sign(w.Exp1) == np.sign(w.Exp2))),
      "| corr %.2f" % np.corrcoef(w.Exp1, w.Exp2)[0, 1])
for e in EXPS:
    d = df[df.exp == e]
    print(e, "median |non-additivity| %.3f vs median |additive| %.3f; synergy(+) %d, antagonism(-) %d" % (
        d.nonadd.abs().median(), d.additive.abs().median(), (d.nonadd > 0).sum(), (d.nonadd < 0).sum()))
