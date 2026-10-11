#!/usr/bin/env python3
"""Step 61: (Supplemental) Table 1 - overlap between reproduced scTARGET-seq targets (main effects, both
experiments) and TARGET targets of the same TFs (Brooks et al. 2019, Supplementary Data 2).
Universe: the 11,980 genes analysed here, minus CHX-responsive genes (Brooks Supplementary Data 1).
Fisher exact test (one-sided, enrichment), BH across TFs; direction agreement on shared targets."""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from common import RES, bh
from xlsx_min import read_xlsx
EXT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "external_Brooks2019")
rows = list(read_xlsx(os.path.join(EXT, "Brooks2019_MOESM5.xlsx")).values())[0][3:]
B = pd.DataFrame([dict(tf=r.get("A"), agi=r.get("B"), target=r.get("C"), lfc=float(r.get("D"))) for r in rows if r.get("C")])
chx = {r.get("A") for r in list(read_xlsx(os.path.join(EXT, "Brooks2019_MOESM4.xlsx")).values())[0][1:] if r.get("A")}
genes = set(np.load(os.path.join(RES, "final", "full_stats.npz"), allow_pickle=True)["genes"])
U = genes - chx
me = pd.read_csv(os.path.join(RES, "consensus_log", "main_effects_all_significant.tsv.gz"), sep="\t")
me = me[me.consensus & me.target.isin(U)]
out = []
for tf in sorted(me.regulator.unique() | set()) if False else sorted(set(pd.read_csv(os.path.join(RES, "consensus_log", "per_TF_summary.tsv"), sep="\t").tf)):
    s = me[me.regulator == tf].set_index("target").beta_meta
    b = B[(B.tf == tf) & B.target.isin(U)].drop_duplicates("target").set_index("target").lfc
    k = len(set(s.index) & set(b.index)); n1, n2, N = len(s), len(b), len(U)
    tab = [[k, n1 - k], [n2 - k, N - n1 - n2 + k]]
    orr, p = stats.fisher_exact(tab, alternative="greater") if n1 and n2 else (np.nan, np.nan)
    sh = list(set(s.index) & set(b.index))
    agree = float(np.mean(np.sign(s[sh]) == np.sign(b[sh]))) if sh else np.nan
    out.append(dict(TF=tf, AGI=B[B.tf == tf].agi.iloc[0], scTARGET_targets=n1, TARGET_targets=n2, shared=k,
                    expected=round(n1 * n2 / N, 1), odds_ratio=orr, p=p,
                    pct_scTARGET_found_by_TARGET=100 * k / n1 if n1 else np.nan,
                    pct_TARGET_recovered_by_scTARGET=100 * k / n2 if n2 else np.nan, same_direction_pct=100 * agree if sh else np.nan))
T = pd.DataFrame(out); T["q_BH"] = bh(T.p.values)
T = T.sort_values("p")
T.to_csv(os.path.join(RES, "paper", "table1_brooks_overlap.tsv"), sep="\t", index=False)
pd.set_option("display.width", 250)
print("universe", len(U), "genes (", len(genes), "analysed minus", len(genes & chx), "CHX-responsive )")
print(T.round({"odds_ratio": 2, "p": 6, "q_BH": 6, "pct_scTARGET_found_by_TARGET": 1, "pct_TARGET_recovered_by_scTARGET": 1, "same_direction_pct": 0}).to_string(index=False))
print("TFs with q<0.05:", int((T.q_BH < 0.05).sum()), "/", len(T), "| p<0.05:", int((T.p < 0.05).sum()))
all_s = set(zip(me.regulator, me.target)); all_b = set(zip(B.tf, B.target))
print("pooled: reproduced edges", len(all_s), "of which in TARGET", len(all_s & all_b))
