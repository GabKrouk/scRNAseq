#!/usr/bin/env python3
"""Step 71: biology of the 461 non-additive interactions (Fig 5 network): GO enrichment and GeneCloud 2.

Gene sets (targets of the 461 TFa x TFb -> gate -> target connections, results/paper/fig5d_candidates.tsv):
  all targets; by gate type (emergent / synergy / competition / reversal); by the direction of the combined
  effect of the pair (both TFs vs neither: induced / repressed); by TF (targets of the gates of each TF);
  hub targets (gated by >= 3 TF pairs).
Universe / background: the 11,980 genes tested in scTARGET-seq (results/final/full_stats.npz).
GO: biological process, TAIR GAF propagated to ancestors (is_a, part_of), hypergeometric, BH over terms with
5-500 annotated universe genes (as step 54). Because scRNA-seq power grows with expression, every term with
q < 0.05 is also tested against 2,000 random gene sets of the same size matched on scTARGET expression
(20 bins of mean log1p CP10k): p_expr_matched.
GeneCloud 2 (G. Krouk; Python rewrite of GeneCloud 2013): words, two-word phrases and UniProt keywords of the
gene annotation (UniProt + NCBI Gene + GO names), hypergeometric vs the same background, BH FDR, redundant
concepts merged. Code and annotation in ~/scdata/genecloud2 (annotation built offline from the TAIR GAF,
go-basic.obo, NCBI gene_info and UniProt UP000006548 files of Oct 2026).
"""
import os, sys, gzip, re, json
import numpy as np, pandas as pd
from scipy import stats, sparse
sys.path.insert(0, os.path.dirname(__file__))
from common import RES, DATA, bh
GC_HOME = os.path.expanduser("~/scdata/genecloud2")
sys.path.insert(0, os.path.join(GC_HOME, "src"))
from genecloud import GeneCloud
from genecloud.compare import draw_compare, compare_table

SRC = "/Users/gabkrouk/Downloads/nrt2q_RNAseq_analysis/genecloud2_sources"
OUT = os.path.join(RES, "go_network"); os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(71)

# ------------------------------------------------------------------ gene sets
P = pd.read_csv(os.path.join(RES, "paper", "fig5d_candidates.tsv"), sep="\t")
P["gate"] = P["shape"].replace({"cancellation": "competition"})
P["a"] = P.pair.str.split("-").str[0]; P["b"] = P.pair.str.split("-").str[1]
U = list(np.load(os.path.join(RES, "final", "full_stats.npz"))["genes"])
sets = {"All network targets": sorted(P.target.unique())}
for g in ("emergent", "synergy", "competition", "reversal"):
    sets[f"Gate: {g}"] = sorted(P.loc[P.gate == g, "target"].unique())
sets["Induced by the pair"] = sorted(P.loc[P.dAB > 0, "target"].unique())
sets["Repressed by the pair"] = sorted(P.loc[P.dAB < 0, "target"].unique())
npair = P.groupby("target").pair.nunique()
sets["Hub targets (>= 3 pairs)"] = sorted(npair[npair >= 3].index)
for t in sorted(set(P.a) | set(P.b)):
    tg = sorted(P.loc[(P.a == t) | (P.b == t), "target"].unique())
    if len(tg) >= 20:
        sets[f"TF: {t}"] = tg
print({k: len(v) for k, v in sets.items()})

# ------------------------------------------------------------------ expression in scTARGET-seq (for matched nulls)
expr = pd.Series(0.0, index=U)
for e in ("Exp1", "Exp2"):
    G = sparse.load_npz(os.path.join(DATA, f"{e}_genes.npz")).tocsr()
    genes = pd.read_csv(os.path.join(DATA, f"{e}_genes.txt"), header=None)[0]
    tot = np.asarray(G.sum(0)).ravel()
    idx = pd.Index(genes).get_indexer(U)
    X = G[idx].multiply(1e4 / tot[None, :]).tocsr(); X.data = np.log1p(X.data)
    expr += np.asarray(X.mean(1)).ravel() / 2
ebin = pd.qcut(expr.rank(method="first"), 20, labels=False)

# ------------------------------------------------------------------ GO (as step 54)
name, ns, par = {}, {}, {}
cur = None
for line in gzip.open(os.path.join(SRC, "go-basic.obo.gz"), "rt"):
    line = line.strip()
    if line == "[Term]": cur = {}; continue
    if line.startswith("["): cur = None; continue
    if cur is None or not line: continue
    if line.startswith("id: GO:"): cur["id"] = line[4:]; par[cur["id"]] = set()
    elif line.startswith("name: ") and "id" in cur: name[cur["id"]] = line[6:]
    elif line.startswith("namespace: ") and "id" in cur: ns[cur["id"]] = line[11:]
    elif line.startswith("is_a: ") and "id" in cur: par[cur["id"]].add(line[6:16])
    elif line.startswith("relationship: part_of ") and "id" in cur: par[cur["id"]].add(line[22:32])
anc = {}
sys.setrecursionlimit(10000)
def ancestors(t):
    if t in anc: return anc[t]
    a = {t}
    for p in par.get(t, ()): a |= ancestors(p)
    anc[t] = a; return a
g2t = {}
for line in gzip.open(os.path.join(SRC, "tair.gaf.gz"), "rt"):
    if line.startswith("!"): continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[8] != "P" or "NOT" in f[3]: continue
    agi = f[1].upper()
    if not re.match(r"AT[1-5CM]G\d{5}$", agi): continue
    g2t.setdefault(agi, set()).update(ancestors(f[4]))
UA = [g for g in U if g in g2t]
cnt_u = {}
for g in UA:
    for t in g2t[g]: cnt_u[t] = cnt_u.get(t, 0) + 1
terms = [t for t, K in cnt_u.items() if ns.get(t) == "biological_process" and 5 <= K <= 500]
tidx = {t: i for i, t in enumerate(terms)}
A = np.zeros((len(UA), len(terms)), bool)
for i, g in enumerate(UA):
    for t in g2t[g]:
        if t in tidx: A[i, tidx[t]] = True
ua_pos = pd.Index(UA); ua_bin = ebin.loc[UA].values


def go_enrich(genes, label):
    S = [g for g in genes if g in ua_pos]
    if len(S) < 5: return pd.DataFrame()
    si = ua_pos.get_indexer(S)
    k = A[si].sum(0); K = A.sum(0); N = len(UA); n = len(S)
    p = stats.hypergeom.sf(k - 1, N, K, n)
    d = pd.DataFrame(dict(set=label, go=terms, term=[name[t] for t in terms], overlap=k, n_set=n, annotated=K,
                          n_universe=N, fold=(k / n) / (K / N), p=p))
    d = d[d.overlap > 0].copy(); d["q_BH"] = bh(d.p.values)
    d["genes"] = [",".join(np.array(S)[A[si, tidx[t]]]) for t in d.go]
    # expression-matched null for q < 0.05 terms
    d["p_expr_matched"] = np.nan
    sig = d.index[d.q_BH < 0.05]
    if len(sig):
        cols = [tidx[t] for t in d.loc[sig, "go"]]
        comp = pd.Series(ua_bin[si]).value_counts()
        pools = {b: np.where(ua_bin == b)[0] for b in comp.index}
        R = 2000; null = np.zeros((R, len(cols)))
        for r in range(R):
            draw = np.concatenate([rng.choice(pools[b], c, replace=False) for b, c in comp.items()])
            null[r] = A[np.ix_(draw, cols)].sum(0)
        d.loc[sig, "p_expr_matched"] = (1 + (null >= d.loc[sig, "overlap"].values[None, :]).sum(0)) / (1 + R)
        d.loc[sig, "expected_expr_matched"] = null.mean(0)
    return d.sort_values("p")


GO = pd.concat([go_enrich(v, k) for k, v in sets.items()], ignore_index=True)
GO.to_csv(os.path.join(OUT, "go_bp_network_sets.tsv"), sep="\t", index=False)
top = GO[(GO.q_BH < 0.05)]
print(top.groupby("set").head(8)[["set", "term", "overlap", "n_set", "fold", "q_BH", "expected_expr_matched", "p_expr_matched"]].round(4).to_string(index=False))

# ------------------------------------------------------------------ GeneCloud 2
gc = GeneCloud(os.path.join(GC_HOME, "arabidopsis.jsonl.gz"), layers=("word", "phrase", "keyword"))
res = {}
for k, v in sets.items():
    r = gc.run(v, background=U, fdr=0.05, min_genes=2)
    res[k] = r
    slug = re.sub(r"[^A-Za-z0-9]+", "_", k).strip("_")
    r.to_tsv(os.path.join(OUT, f"genecloud_{slug}.tsv"))
    gc.plot(r, os.path.join(OUT, f"genecloud_{slug}.pdf"), title=f"{k} ({len(r.study)} genes)")
    gc.plot(r, os.path.join(OUT, f"genecloud_{slug}.png"), title=f"{k} ({len(r.study)} genes)")
    gc.html(r, os.path.join(OUT, f"genecloud_{slug}.html"), title=f"{k} ({len(r.study)} genes)")
    s = r.significant
    print(f"{k}: {len(r.study)} genes, {len(s)} significant concepts:",
          ", ".join(f"{a} ({b}/{c}, FDR {d:.1g})" for a, b, c, d in zip(s.display.head(12), s.k.head(12), s.n.head(12), s.fdr.head(12))))
main = {k: res[k] for k in sets if not k.startswith("TF:")}
draw_compare(main, os.path.join(OUT, "genecloud_compare_sets.pdf"), per_set=6, title="GeneCloud 2 - targets of the 461 interactions")
draw_compare({k: res[k] for k in sets if k.startswith("TF:")}, os.path.join(OUT, "genecloud_compare_TFs.pdf"), per_set=5,
             title="GeneCloud 2 - targets of the gates of each TF")
compare_table(res, per_set=8).to_csv(os.path.join(OUT, "genecloud_compare_table.tsv"), sep="\t", index=False)
json.dump({k: len(v) for k, v in sets.items()}, open(os.path.join(OUT, "set_sizes.json"), "w"), indent=1)
