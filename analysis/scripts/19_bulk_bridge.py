#!/usr/bin/env python3
"""Step 19: bridge single-cell (Exp1, Exp2) and the ERF5 x CRF4 bulk DoubleTARGET
RNA-seq (Fig 7; Will's counts.txt, 2 experiments x EV/ERF5/CRF4/ERF5+CRF4 x 3).
Bulk model per gene (log2 CPM, E2_EC.2 and E2_E.1 removed as in Will's Rmd):
    y ~ experiment + ERF5 + CRF4 + ERF5:CRF4   (OLS)
Compared with the single-cell ERF5 / CRF4 main effects and the ERF5 x CRF4
interaction (step 5), and Will's synergy/additive gene lists vs the
single-cell combinatorial targets (step 13/14)."""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from common import *
from xlsx_min import read_xlsx

W = "/Users/gabkrouk/Library/CloudStorage/Dropbox/Science/Papers/Emanuel New TARGET paper/Figure Will"
OUT = os.path.join(RES, "bulk_bridge"); os.makedirs(OUT, exist_ok=True)
C = pd.read_csv(os.path.join(W, "counts.txt"), sep=" ", index_col=0)
C = C.drop(columns=["E2_EC.2", "E2_E.1"])
cpm = C / C.sum() * 1e6
keep = (cpm > 1).sum(1) >= 3
L = np.log2(cpm[keep] + 1)
s = pd.Series(L.columns)
exp2 = s.str.startswith("E2").astype(float).values
con = s.str.extract(r"_(\w+)\.")[0]
E = con.isin(["E", "EC"]).astype(float).values; Cc = con.isin(["C", "EC"]).astype(float).values
D = np.column_stack([np.ones(len(s)), exp2, E, Cc, E * Cc])
Y = L.values.T
B = np.linalg.lstsq(D, Y, rcond=None)[0]; R = Y - D @ B
dfr = len(s) - D.shape[1]; s2 = (R ** 2).sum(0) / dfr
cv = np.linalg.inv(D.T @ D)
bulk = pd.DataFrame(index=L.index)
for k, nm in ((2, "ERF5"), (3, "CRF4"), (4, "inter")):
    t = B[k] / np.sqrt(s2 * cv[k, k])
    bulk[f"bulk_{nm}"] = B[k]; bulk[f"bulk_p_{nm}"] = 2 * stats.t.sf(np.abs(t), dfr)
bulk.to_csv(os.path.join(OUT, "bulk_ERF5xCRF4_effects.tsv"), sep="\t")
print("bulk genes:", len(bulk), "| interaction p<0.05:", (bulk.bulk_p_inter < .05).sum(),
      "| ERF5 p<0.05:", (bulk.bulk_p_ERF5 < .05).sum(), "| CRF4 p<0.05:", (bulk.bulk_p_CRF4 < .05).sum())

st = np.load(os.path.join(RES, "consensus_ge1", "full_stats.npz"))
genes = list(st["genes"]); tfs = list(st["tfs"]); pairs = list(st["pairs"])
sc = pd.DataFrame(index=genes)
for e in EXPS:
    for tf in ("ERF5", "CRF4"):
        sc[f"{e}_{tf}"] = st[f"{e}_main_beta"][tfs.index(tf)]
        sc[f"{e}_se_{tf}"] = st[f"{e}_main_se"][tfs.index(tf)]
    k = pairs.index("CRF4-ERF5")
    sc[f"{e}_inter"] = st[f"{e}_inter_beta"][k]; sc[f"{e}_se_inter"] = st[f"{e}_inter_se"][k]
for nm in ("ERF5", "CRF4", "inter"):
    w1, w2 = 1 / sc[f"Exp1_se_{nm}"] ** 2, 1 / sc[f"Exp2_se_{nm}"] ** 2
    sc[f"meta_{nm}"] = (sc[f"Exp1_{nm}"] * w1 + sc[f"Exp2_{nm}"] * w2) / (w1 + w2)
M = sc.join(bulk, how="inner")
print("genes in both:", len(M))
rows = []
for nm in ("ERF5", "CRF4", "inter"):
    sel = M[f"bulk_p_{nm}"] < 0.05
    for src in ("Exp1", "Exp2", "meta"):
        r_all = stats.spearmanr(M[f"{src}_{nm}"], M[f"bulk_{nm}"]).correlation
        r_sel = stats.spearmanr(M.loc[sel, f"{src}_{nm}"], M.loc[sel, f"bulk_{nm}"]).correlation
        agree = np.mean(np.sign(M.loc[sel, f"{src}_{nm}"]) == np.sign(M.loc[sel, f"bulk_{nm}"]))
        bt = stats.binomtest(int(agree * sel.sum()), int(sel.sum()), 0.5).pvalue
        rows.append(dict(effect=nm, singlecell=src, rho_all_genes=r_all, n_bulk_p05=int(sel.sum()),
                         rho_bulk_DEGs=r_sel, sign_agreement_bulk_DEGs=agree, binom_p=bt))
cmp_ = pd.DataFrame(rows); cmp_.to_csv(os.path.join(OUT, "sc_vs_bulk_effects.tsv"), sep="\t", index=False)
pd.set_option("display.width", 200); print(cmp_.round(4).to_string(index=False))

# Will's gene lists vs single-cell combinatorial targets
def genes_of(f):
    rows = read_xlsx(os.path.join(W, f))["Sheet 1"]
    return {r.get("A") for r in rows[1:] if r.get("A")}
lists = {"Will synergy (ERF5xCRF4 interaction p<0.05)": genes_of("Synergy_Consensus_Clusters.xlsx"),
         "Will additive DEGs": genes_of("Additive_Clusters.xlsx")}
comb = pd.read_csv(os.path.join(RES, "combinatorial", "four_state_effects_all_consensus.tsv"), sep="\t")
prop = pd.read_csv(os.path.join(RES, "combinatorial", "proper_combinatorial_targets.tsv"), sep="\t")
cons_main = pd.read_csv(os.path.join(RES, "consensus_ge1", "main_effects_all_significant.tsv.gz"), sep="\t")
cons_main = cons_main[cons_main.consensus]
U = set(M.index)
scsets = {"sc consensus interaction, any pair": set(comb.target),
          "sc proper combinatorial, any pair": set(prop.target),
          "sc proper combinatorial, ERF5 pairs": set(prop[(prop.TFa == "ERF5") | (prop.TFb == "ERF5")].target),
          "sc consensus TF->gene, any TF": set(cons_main.target),
          "sc consensus ERF5 or CRF4 target": set(cons_main[cons_main.regulator.isin(["ERF5", "CRF4"])].target)}
ov = []
for ln, lg in lists.items():
    for sn, sg in scsets.items():
        a = len(lg & sg & U); b = len((lg & U) - sg); c = len((sg & U) - lg); d = len(U) - a - b - c
        orr, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
        ov.append(dict(bulk_list=ln, singlecell_set=sn, n_bulk=len(lg & U), n_sc=len(sg & U), overlap=a,
                       expected=round(len(lg & U) * len(sg & U) / len(U), 1), odds_ratio=orr, p=p))
ov = pd.DataFrame(ov); ov.to_csv(os.path.join(OUT, "overlap_bulk_lists_vs_singlecell.tsv"), sep="\t", index=False)
print(ov.round(4).to_string(index=False))
M.to_csv(os.path.join(OUT, "sc_and_bulk_per_gene.tsv"), sep="\t")
