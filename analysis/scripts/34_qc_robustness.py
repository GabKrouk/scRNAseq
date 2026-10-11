#!/usr/bin/env python3
"""Step 34: (a) TF-detection sensitivity vs sequencing depth from raw Exp2 counts
(Reviewer 2 #2, Reviewer 3 #3); (b) robustness of the consensus network to strict
cell QC: re-run step 5 (main effects + all TF x TF interactions, both experiments)
keeping only cells that pass raw-count QC (no low-UMI / low-gene MAD flag, mito <5%,
plastid <5%, not a called doublet) and that are present in the CellRanger output."""
import os, sys, json
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
os.environ["N_PERM"] = "0"
src = open(os.path.join(os.path.dirname(__file__), "05_consensus.py")).read()
src = src[:src.index("t0 = time.time()")]
exec(compile(src, "05", "exec"))
OUT = os.path.join(RES, "qc_raw")

# ---------- (a) TF detection vs depth, Exp2 raw
q2 = pd.read_csv(os.path.join(OUT, "Exp2_cell_qc.tsv.gz"), sep="\t")
q2["log10_umi"] = np.log10(q2.umi)
bins = pd.qcut(q2.log10_umi, 10)
det = q2.groupby(bins, observed=True).agg(umi_median=("umi", "median"), frac_no_TF=("tf_reads", lambda x: np.mean(x == 0)),
                                          median_n_TF=("n_tf", "median"), frac_RFP=("RFP", lambda x: np.mean(x > 0)), n=("umi", "size")).reset_index(drop=True)
det.to_csv(os.path.join(OUT, "Exp2_TF_detection_vs_depth.tsv"), sep="\t", index=False)
print(det.round(3).to_string(index=False))
import statsmodels.api as sm
m = sm.GLM((q2.tf_reads == 0).astype(float), sm.add_constant(q2.log10_umi), family=sm.families.Binomial()).fit()
print("logistic P(no TF read) ~ log10 UMI: slope %.2f (p=%.1e)" % (m.params.iloc[1], m.pvalues.iloc[1]))
for u in (1e4, 2e4, 5e4, 1e5):
    print("  predicted fraction without TF read at %d UMI: %.3f" % (u, m.predict(np.array([[1, np.log10(u)]]))[0]))

# ---------- (b) strict-QC re-run of the consensus analysis
keep = {}
for e in EXPS:
    q = pd.read_csv(os.path.join(OUT, f"{e}_cell_qc.tsv.gz"), sep="\t")
    bad = q.low_umi | q.low_genes | (q.pct_mito > 5) | (q.pct_plastid > 5) | q.doublet
    good = set(q.barcode[~bad])
    cells = ds[e]["cells"].cell.str.replace(".", "-", regex=False)
    keep[e] = cells.isin(good).values
    print(e, "cells in analysis:", len(cells), "-> strict QC:", int(keep[e].sum()))
F0 = {e: fit_exp(X[e], Z[e], Y[e]) for e in EXPS}
F1 = {e: fit_exp(X[e][keep[e]], Z[e][keep[e]][:, Z[e][keep[e]].std(0) > 0], Y[e][keep[e]]) for e in EXPS}
s0, s1 = summarise(F0, "all"), summarise(F1, "strict")
res = {}
for kind in ("main", "inter"):
    m0, m1 = meta(F0["Exp1"][kind], F0["Exp2"][kind]), meta(F1["Exp1"][kind], F1["Exp2"][kind])
    c0, c1 = consensus(F0["Exp1"][kind], F0["Exp2"][kind], m0), consensus(F1["Exp1"][kind], F1["Exp2"][kind], m1)
    ok = ~np.isnan(m0["beta"]) & ~np.isnan(m1["beta"])
    res[kind] = dict(consensus_all_cells=int(c0.sum()), consensus_strict_qc=int(c1.sum()), retained=int((c0 & c1).sum()),
                     retained_frac=float((c0 & c1).sum() / max(c0.sum(), 1)),
                     corr_meta_beta=float(np.corrcoef(m0["beta"][ok], m1["beta"][ok])[0, 1]),
                     same_sign_of_original_consensus=float(np.mean(np.sign(m0["beta"][c0]) == np.sign(m1["beta"][c0]))))
print(json.dumps(res, indent=1))
json.dump(res, open(os.path.join(OUT, "strict_qc_robustness.json"), "w"), indent=1)
