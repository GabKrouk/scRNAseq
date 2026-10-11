#!/usr/bin/env python3
"""Step 57 (R1 M5): wounding/protoplasting stress markers in cells with no detected GR-TF vs cells
with GR-TFs, and vs the number of GR-TFs, per experiment (QC cells, log1p CP10k)."""
import os, sys, json
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from common import *
MK = {"ZAT12": "AT5G59820", "WRKY33": "AT2G38470", "ZAT10": "AT1G27730", "ERF104": "AT5G61600", "HSP17.6A": "AT5G12030",
      "HSP70": "AT3G12580", "CYP81D8": "AT4G37370", "WRKY40": "AT1G80840", "MYB51": "AT1G18570", "JAZ1": "AT1G19180"}
ds = {e: load(e) for e in EXPS}; _, tfs = analysis_sets(ds)
rows = []
for e in EXPS:
    d = ds[e]; n = (d["T"][tfs].values >= 1).sum(1)
    for nm, g in MK.items():
        if g not in d["genes"]: continue
        y = response(d, [g])[:, 0]
        rho = stats.spearmanr(n, y).correlation
        rows.append(dict(exp=e, marker=nm, agi=g, frac_expr=float((y > 0).mean()), mean_noTF=float(y[n == 0].mean()),
                         mean_TF=float(y[n > 0].mean()), p_MW=float(stats.mannwhitneyu(y[n == 0], y[n > 0]).pvalue), rho_nTF=float(rho)))
df = pd.DataFrame(rows); print(df.round(3).to_string(index=False))
df.to_csv(os.path.join(RES, "paper", "stress_markers.tsv"), sep="\t", index=False)
