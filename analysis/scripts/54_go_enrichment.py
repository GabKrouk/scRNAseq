#!/usr/bin/env python3
"""Step 54: GO biological-process enrichment (TAIR GAF + go-basic.obo, annotations propagated to
ancestors via is_a/part_of; hypergeometric test; BH across terms with 5-500 annotated genes; universe =
expressed genes of the dataset). Gene sets: bulk ERF5 x CRF4 classes (Fig 6) and the single-cell
non-additive interaction targets (Fig 5)."""
import os, sys, gzip, re
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from common import RES, bh
SRC = "/Users/gabkrouk/Downloads/nrt2q_RNAseq_analysis/genecloud2_sources"
OUT = os.path.join(RES, "go"); os.makedirs(OUT, exist_ok=True)
name, ns, par = {}, {}, {}
cur = None
for line in gzip.open(os.path.join(SRC, "go-basic.obo.gz"), "rt"):
    line = line.strip()
    if line == "[Term]": cur = {}; continue
    if line.startswith("[") : cur = None; continue
    if cur is None or not line: continue
    if line.startswith("id: GO:"): cur["id"] = line[4:]; par[cur["id"]] = set()
    elif line.startswith("name: ") and "id" in cur: name[cur["id"]] = line[6:]
    elif line.startswith("namespace: ") and "id" in cur: ns[cur["id"]] = line[11:]
    elif line.startswith("is_a: ") and "id" in cur: par[cur["id"]].add(line[6:16])
    elif line.startswith("relationship: part_of ") and "id" in cur: par[cur["id"]].add(line[22:32])
anc = {}
def ancestors(t):
    if t in anc: return anc[t]
    a = {t}
    for p in par.get(t, ()): a |= ancestors(p)
    anc[t] = a; return a
sys.setrecursionlimit(10000)
g2t = {}
for line in gzip.open(os.path.join(SRC, "tair.gaf.gz"), "rt"):
    if line.startswith("!"): continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[8] != "P" or "NOT" in f[3]: continue
    agi = f[1].upper()
    if not re.match(r"AT[1-5CM]G\d{5}$", agi): continue
    g2t.setdefault(agi, set()).update(ancestors(f[4]))
def enrich(genes, universe, label):
    U = [g for g in universe if g in g2t]; S = set(genes) & set(U)
    if len(S) < 5: return pd.DataFrame()
    cnt_u, cnt_s = {}, {}
    for g in U:
        for t in g2t[g]: cnt_u[t] = cnt_u.get(t, 0) + 1
    for g in S:
        for t in g2t[g]: cnt_s[t] = cnt_s.get(t, 0) + 1
    rows = []
    for t, k in cnt_s.items():
        K = cnt_u[t]
        if ns.get(t) != "biological_process" or K < 5 or K > 500: continue
        rows.append(dict(set=label, go=t, term=name.get(t, t), n_set=len(S), n_universe=len(U), overlap=k, annotated=K,
                         fold=(k / len(S)) / (K / len(U)), p=stats.hypergeom.sf(k - 1, len(U), K, len(S))))
    d = pd.DataFrame(rows)
    if len(d): d["q_BH"] = bh(d.p.values); d = d.sort_values("p")
    return d
res = []
cl = pd.read_csv(os.path.join(RES, "bulk_bridge", "Fig6_gene_classes.tsv"), sep="\t")
allb = pd.read_csv(os.path.join(RES, "bulk_bridge", "bulk_condition_contrasts_per_gene.tsv"), sep="\t", index_col=0)
Ub = list(allb.index)
for c in cl["class"].unique():
    sub = cl[cl["class"] == c]
    up = sub.log2FC_CRF4 if c != "combination-specific" else sub.log2FC_ERF5_CRF4
    for d, m in (("up", up > 0), ("down", up < 0)):
        res.append(enrich(sub.gene[m], Ub, f"bulk {c} {d}"))
sl = pd.read_csv(os.path.join(RES, "combinatorial_push", "supinfo_461_list.tsv"), sep="\t")
Us = [l.strip() for l in open(os.path.join(RES, "..", "data", "Exp1_genes.txt"))]
me = pd.read_csv(os.path.join(RES, "consensus_log", "main_effects_all_significant.tsv.gz"), sep="\t")
import numpy as _np
Us = list(_np.load(os.path.join(RES, "final", "full_stats.npz"), allow_pickle=True)["genes"])
res.append(enrich(sl.target.unique(), Us, "single-cell non-additive targets (461)"))
for tf in ("NAC4", "CRF4"):
    m = me[(me.regulator == tf) & me.consensus]
    for d, s in (("up", m.beta_meta > 0), ("down", m.beta_meta < 0)):
        res.append(enrich(m.target[s], Us, f"single-cell {tf} {d}"))
R = pd.concat([r for r in res if len(r)])
R.to_csv(os.path.join(OUT, "go_enrichment_all.tsv"), sep="\t", index=False)
top = R[R.q_BH < 0.05].groupby("set").head(6)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
print(R.groupby("set").apply(lambda x: pd.Series(dict(n_set=x.n_set.iloc[0], n_q05=(x.q_BH < .05).sum()))))
print(top[["set", "term", "overlap", "fold", "q_BH"]].to_string(index=False))
