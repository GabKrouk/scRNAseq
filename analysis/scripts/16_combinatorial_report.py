#!/usr/bin/env python3
"""Step 16: self-contained HTML for the combinatorial (4-state) analysis."""
import os, base64
import pandas as pd
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FG = os.path.join(B, "figures", "combinatorial"); D = os.path.join(B, "results", "combinatorial")
img = lambda n, w="100%": f'<img src="data:image/png;base64,{base64.b64encode(open(os.path.join(FG, n), "rb").read()).decode()}" style="max-width:{w}">'
P = pd.read_csv(os.path.join(D, "proper_combinatorial_targets.tsv"), sep="\t")
core = P[P.core]
top = core.head(25)
rows = "".join(f"<tr><td>{r.TFa} &times; {r.TFb}</td><td>{r.target}</td><td>{r.shape}</td><td>{r.a_meta:+.2f}</td><td>{r.b_meta:+.2f}</td><td>{r.both_meta:+.2f}</td><td>{r.p_interaction_meta:.1e}</td></tr>" for r in top.itertuples())
html = f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>TF pair combinatorial targets</title>
<style>body{{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;max-width:1100px;margin:0 auto;padding:24px 16px 60px;color:#1d1d1f;background:#fff;line-height:1.5}}
h2{{border-bottom:2px solid #e3e3e3;padding-bottom:.2em;margin-top:2em}} table{{border-collapse:collapse;font-size:.88em;width:100%}} th,td{{border-bottom:1px solid #e3e3e3;padding:4px 8px;text-align:left}} th{{background:#f6f7f8}}
.box{{background:#f6f7f8;border-left:4px solid #CC3311;padding:10px 16px}} code{{background:#f6f7f8;padding:1px 4px}}</style></head><body>
<h1>TFa &times; TFb &rarr; target: combinatorial effects that the single TFs cannot predict</h1>
<div class="box"><b>Funnel.</b> 1,159 reproducible TF&times;TF interactions &rarr; <b>366 proper combinatorial cases</b> ({P.target.nunique()} genes, {P.pair.nunique()} pairs) &rarr; <b>{len(core)} core cases</b> ({core.target.nunique()} genes, {core.pair.nunique()} pairs).
<br><b>Proper</b> means four things. First, the "both" state differs from "a alone" and from "b alone" at p&lt;0.05 <i>in each experiment</i>, in the same direction. Second, the interaction term is larger than the larger single-TF effect. Third, "both" and the interaction have the same sign in Exp1 and Exp2. Fourth, the single-TF effects do not contradict each other between experiments.
<br><b>Core</b> means two more conditions. The case is pair-specific: on that gene, neither TF is combinatorial with more than 2 partners. It also stays proper when the number of other TFs detected in the cell is adjusted as categories (one dummy per count).</div>
<p>Every state effect is adjusted for depth, RFP, GR-TF load (non-linear), library, cell type and the other 21 TFs. The bars pool both experiments; the circles show Exp1 and Exp2 separately.</p>
<h2>Overview</h2>{img("overview_shapes.png","48%")} {img("pair_matrix.png","50%")}{img("overview_profiles.png")}
<p><b>AND-gate:</b> neither TF alone changes the gene, but the two together do. <b>Reversal:</b> a single TF pushes the gene one way and the pair pushes it the other way. <b>Cancellation:</b> the single TFs act, but together they return to the "neither" level. <b>Synergy:</b> both singles push the same way and the pair goes far beyond. Almost no case is pure same-direction synergy: the reproducible combinatorial logic is mostly AND-gates and reversals.</p>
<h2>Best core examples</h2>{img("gallery_best_by_shape.png")}
<h2>Top 25 core cases (meta-analysis effects vs "neither", log expression)</h2>
<table><tr><th>Pair</th><th>Target</th><th>Shape</th><th>a alone</th><th>b alone</th><th>both</th><th>interaction p</th></tr>{rows}</table>
<h2>Files</h2><p>All 366 cases, 12 per page, Exp1 / Exp2 / pooled: <code>figures/combinatorial/ALL_proper_combinatorial_targets.pdf</code>. Tables: <code>results/combinatorial/proper_combinatorial_targets.tsv</code> (core and pair-specific flags), <code>four_state_effects_all_consensus.tsv</code> (all 1,159), <code>state_means_for_plots.tsv</code>. Scripts 13&ndash;16.</p></body></html>"""
open(os.path.join(B, "REPORT_combinatorial_targets.html"), "w").write(html); print("ok")
