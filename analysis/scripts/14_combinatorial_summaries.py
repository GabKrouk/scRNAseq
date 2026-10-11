#!/usr/bin/env python3
"""Step 14: per-state summaries for plotting + pair specificity.
For each (pair, target): mean adjusted expression and 95% CI per state
(neither / TFa only / TFb only / both), expressed as change from 'neither',
in Exp1, Exp2 and both combined (cells pooled after centring each experiment
on its own 'neither' mean).  Pair specificity = number of different partners
with which the same TF shows a proper combinatorial effect on the same target."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import RES
OUT = os.path.join(RES, "combinatorial")
df = pd.read_csv(os.path.join(OUT, "four_state_effects_all_consensus.tsv"), sep="\t")
cells = pd.read_csv(os.path.join(OUT, "cells_adjusted_expression.tsv.gz"), sep="\t")
ref = cells[cells.state == "neither"].groupby(["pair", "target", "exp"]).adj_expr.mean().rename("ref")
cells = cells.join(ref, on=["pair", "target", "exp"])
cells["delta"] = cells.adj_expr - cells.ref
both = cells.copy(); both["exp"] = "Combined"
allc = pd.concat([cells, both])
g = allc.groupby(["pair", "target", "exp", "state"]).delta
s = pd.DataFrame({"mean": g.mean(), "sd": g.std(), "n": g.size()}).reset_index()
s["ci"] = 1.96 * s["sd"] / np.sqrt(s["n"])
s.to_csv(os.path.join(OUT, "state_means_for_plots.tsv"), sep="\t", index=False)

prop = df[df.proper_combinatorial].copy()
# robustness: still 'proper' when the number of other TFs detected is adjusted as categories
rob = pd.read_csv(os.path.join(RES, "combinatorial_nTFcat", "four_state_effects_all_consensus.tsv"), sep="\t")
rob = set(zip(rob[rob.proper_combinatorial].pair, rob[rob.proper_combinatorial].target))
prop["robust_to_nTF"] = [(p_, t_) in rob for p_, t_ in zip(prop.pair, prop.target)]
long = pd.concat([prop.assign(tf=prop.TFa, partner=prop.TFb), prop.assign(tf=prop.TFb, partner=prop.TFa)])
npart = long.groupby(["target", "tf"]).partner.nunique().rename("n_partners").reset_index()
spec = prop.merge(npart.rename(columns={"tf": "TFa", "n_partners": "partners_TFa"}), on=["target", "TFa"]) \
           .merge(npart.rename(columns={"tf": "TFb", "n_partners": "partners_TFb"}), on=["target", "TFb"])
spec["pair_specific"] = (spec.partners_TFa <= 2) & (spec.partners_TFb <= 2)
spec["core"] = spec.pair_specific & spec.robust_to_nTF
spec = spec.sort_values(["core", "score"], ascending=[False, False])
spec.to_csv(os.path.join(OUT, "proper_combinatorial_targets.tsv"), sep="\t", index=False)
print("proper:", len(spec), "targets", spec.target.nunique())
print("pair-specific (each TF combinatorial with <=2 partners on that target):", spec.pair_specific.sum(),
      "targets", spec[spec.pair_specific].target.nunique(), "pairs", spec[spec.pair_specific].pair.nunique())
print("robust to nTF categories:", spec.robust_to_nTF.sum(), "| core (robust & pair-specific):", spec.core.sum(), "genes", spec[spec.core].target.nunique(), "pairs", spec[spec.core].pair.nunique())
print(pd.crosstab(spec["shape"], [spec.pair_specific, spec.robust_to_nTF]))
print(spec[spec.core].pair.value_counts().head(12).to_string())
print(spec[spec.pair_specific].pair.value_counts().head(15).to_string())
