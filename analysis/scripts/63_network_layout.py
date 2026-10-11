#!/usr/bin/env python3
"""Step 63: network of the 461 reproduced non-additive interactions for Fig 5.
Nodes: TFs; 'gates' = one per TF pair and interaction shape (emergent / synergy / competition / reversal);
target genes. Edges: TFa -> gate, TFb -> gate, gate -> target. Layout: Fruchterman-Reingold (numpy, fixed seed)."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import RES
OUT = os.path.join(RES, "paper")
P = pd.read_csv(os.path.join(OUT, "fig5d_candidates.tsv"), sep="\t")
lab = {"emergent": "emergent", "synergy": "synergy", "cancellation": "competition", "reversal": "reversal"}
P["kind"] = P["shape"].map(lab)
P["a"] = P.pair.str.split("-").str[0]; P["b"] = P.pair.str.split("-").str[1]
P["gate"] = P.pair + "|" + P.kind
nodes = {}
def add(n, t):
    if n not in nodes: nodes[n] = t
for _, r in P.iterrows():
    add("TF:" + r.a, "TF"); add("TF:" + r.b, "TF"); add("G:" + r.gate, r.kind); add("T:" + r.target, "target")
edges = set()
for _, r in P.iterrows():
    edges.add(("TF:" + r.a, "G:" + r.gate)); edges.add(("TF:" + r.b, "G:" + r.gate)); edges.add(("G:" + r.gate, "T:" + r.target))
names = list(nodes); idx = {n: i for i, n in enumerate(names)}
E = np.array([(idx[a], idx[b]) for a, b in edges])
n = len(names); rng = np.random.default_rng(63)
pos = rng.normal(size=(n, 2))
k = 1.0 / np.sqrt(n) * 2.2
w = np.array([3.0 if names[i].startswith("TF:") else 1.0 for i in range(n)])     # heavier TF nodes repel more
t = 0.15
for it in range(600):
    d = pos[:, None, :] - pos[None, :, :]
    dist = np.sqrt((d ** 2).sum(-1)) + 1e-9
    rep = (k * k / dist ** 2)[:, :, None] * d * (w[:, None] * w[None, :])[:, :, None]
    disp = rep.sum(1)
    dv = pos[E[:, 0]] - pos[E[:, 1]]; dl = np.sqrt((dv ** 2).sum(1)) + 1e-9
    att = (dl / k)[:, None] * dv
    np.add.at(disp, E[:, 0], -att); np.add.at(disp, E[:, 1], att)
    disp -= 0.02 * pos                                                               # weak gravity
    L = np.sqrt((disp ** 2).sum(1)) + 1e-9
    pos += disp / L[:, None] * np.minimum(L, t)[:, None]
    t = max(0.005, t * 0.992)
pos = (pos - pos.mean(0)) / pos.std()
# second pass: TFs fixed on a circle (order taken from the free layout), gates and targets inside it
isTF = np.array([x.startswith("TF:") for x in names])
ang = np.arctan2(pos[isTF, 1], pos[isTF, 0]); order = np.argsort(ang)
tfi = np.where(isTF)[0][order]; theta = np.linspace(0, 2 * np.pi, len(tfi), endpoint=False)
RAD = 2.0
fixed = np.zeros((n, 2)); fixed[tfi] = RAD * np.column_stack([np.cos(theta), np.sin(theta)])
pos[isTF] = fixed[isTF]; pos[~isTF] *= 0.8
t = 0.1
for it in range(500):
    d = pos[:, None, :] - pos[None, :, :]
    dist = np.sqrt((d ** 2).sum(-1)) + 1e-9
    disp = ((k * k / dist ** 2)[:, :, None] * d).sum(1)
    dv = pos[E[:, 0]] - pos[E[:, 1]]; dl = np.sqrt((dv ** 2).sum(1)) + 1e-9
    att = (dl / k)[:, None] * dv
    np.add.at(disp, E[:, 0], -att); np.add.at(disp, E[:, 1], att)
    L = np.sqrt((disp ** 2).sum(1)) + 1e-9
    pos += disp / L[:, None] * np.minimum(L, t)[:, None]
    pos[isTF] = fixed[isTF]
    r = np.sqrt((pos[~isTF] ** 2).sum(1)); over = r > 0.88 * RAD
    pos[np.where(~isTF)[0][over]] *= (0.88 * RAD / r[over])[:, None]
    t = max(0.004, t * 0.99)
deg = pd.Series(0, index=names, dtype=int)
for a, b in edges: deg[a] += 1; deg[b] += 1
gsize = P.groupby("gate").size()
N = pd.DataFrame({"node": names, "type": [nodes[x] for x in names], "x": pos[:, 0], "y": pos[:, 1], "degree": deg.values})
N["label"] = N.node.str.replace(r"^(TF|G|T):", "", regex=True)
N["n_targets"] = [gsize.get(l, np.nan) if t not in ("TF", "target") else np.nan for l, t in zip(N.label, N.type)]
N.to_csv(os.path.join(OUT, "network_nodes.tsv"), sep="\t", index=False)
pd.DataFrame([(a, b) for a, b in edges], columns=["from", "to"]).to_csv(os.path.join(OUT, "network_edges.tsv"), sep="\t", index=False)
print(N.type.value_counts().to_dict(), "edges", len(edges))
print(N[N.type == "TF"].sort_values("degree", ascending=False)[["label", "degree"]].head(8).to_string(index=False))
