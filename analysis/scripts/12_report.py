#!/usr/bin/env python3
"""Step 12: assemble the self-contained HTML report (figures embedded)."""
import os, base64, json
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
F = lambda n: os.path.join(BASE, "figures", n + ".png")
R = lambda *p: os.path.join(BASE, "results", *p)


def img(name, width="100%"):
    b = base64.b64encode(open(F(name), "rb").read()).decode()
    return f'<img src="data:image/png;base64,{b}" style="max-width:{width}">'


pt = pd.read_csv(R("consensus_ge1", "per_TF_summary.tsv"), sep="\t").sort_values("n_consensus", ascending=False)
pp = pd.read_csv(R("consensus_ge1", "per_pair_summary.tsv"), sep="\t")
sm = json.load(open(R("consensus_ge1", "summary.json")))
tm = pd.read_csv(R("transfer_model_comparison.tsv"), sep="\t", skiprows=3, header=None,
                 names=["model", "r12", "r21", "n12", "n21"])
dose = pd.read_csv(R("consensus_ge1", "interaction_dose_sensitivity.tsv"), sep="\t")

tf_rows = "".join(
    f"<tr><td>{r.tf}</td><td>{r.cells_pos_Exp1}</td><td>{r.cells_pos_Exp2}</td><td>{r.n_q05_Exp1}</td>"
    f"<td>{r.n_q05_Exp2}</td><td>{r.r_t_Exp1_Exp2:.2f}</td><td><b>{r.n_consensus}</b></td>"
    f"<td>{r.n_consensus_up}/{r.n_consensus_down}</td></tr>" for r in pt.itertuples())
top_pairs = pp.sort_values("n_consensus_interaction", ascending=False).head(12)
pair_rows = "".join(
    f"<tr><td>{r['pair']}</td><td>{r['n_++_Exp1']}/{r['n_++_Exp2']}</td><td>{r['n_inter_q05_Exp1']}</td><td>{r['n_inter_q05_Exp2']}</td>"
    f"<td><b>{r['n_consensus_interaction']}</b></td><td>{r['n_consensus_cotargets']}</td></tr>"
    for _, r in top_pairs.iterrows())
paper_pairs = pp[pp.pair.isin(["CRF4-ERF5", "ERF5-NAC4", "ERF5-LBD37", "ERF5-HHO2", "bZIP3-ERF5", "ERF5-HB6", "bZIP3-NAC4", "ERF5-WRKY54"])]
paper_rows = "".join(
    f"<tr><td>{r['pair']}</td><td>{r['n_++_Exp1']}/{r['n_++_Exp2']}</td><td>{r['n_consensus_interaction']}</td><td>{r['n_consensus_cotargets']}</td></tr>"
    for _, r in paper_pairs.sort_values("n_consensus_interaction", ascending=False).iterrows())

html = f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>scTARGET-seq reproducibility</title>
<style>
:root{{--fg:#1d1d1f;--muted:#5f6368;--bg:#fff;--line:#e3e3e3;--accent:#D55E00;--soft:#f6f7f8}}
body{{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;color:var(--fg);background:var(--bg);max-width:980px;margin:0 auto;padding:24px 16px 80px;line-height:1.5}}
h1{{font-size:1.7em;margin-bottom:.2em}} h2{{margin-top:2.2em;border-bottom:2px solid var(--line);padding-bottom:.2em}} h3{{margin-top:1.6em}}
.sub{{color:var(--muted)}} table{{border-collapse:collapse;font-size:.88em;margin:.6em 0;width:100%}} th,td{{border-bottom:1px solid var(--line);padding:4px 8px;text-align:left}}
th{{background:var(--soft)}} .box{{background:var(--soft);border-left:4px solid var(--accent);padding:10px 16px;margin:1em 0}}
.warn{{border-left-color:#0072B2}} code{{background:var(--soft);padding:1px 4px;border-radius:3px}} figure{{margin:1.2em 0}} figcaption{{color:var(--muted);font-size:.9em}}
.tw{{overflow-x:auto}}
</style></head><body>
<h1>scTARGET-seq: combining Exp1 and Exp2, and what reproduces</h1>
<p class="sub">Revision analyses for TPC-2026-0867 (Emanuel, Carluer, Hinckley et al.). Generated 2026-10-04 on gab's Mac from
<code>norm_matrix-Experiment1.txt</code> (2,313 cells) and <code>norm_counts_reseq_Experiment2.txt</code> (6,374 cells, 6 libraries).
Input data were not modified. Every number below is produced by the scripts in <code>reproducibility/scripts/</code>.</p>

<div class="box"><b>Bottom line</b><ul>
<li><b>TF identity reproducibly changes the transcriptome.</b> After adjusting for depth, GR-TF load, RFP, library and cell identity, <b>1,972 TF&rarr;gene effects</b> (on 1,696 genes) are significant in the combined analysis and reproduce, with the same sign, in each experiment. In 50 label permutations the same procedure finds 0.24 such edges on average (max 2), an empirical FDR of 0.01%.</li>
<li><b>Models trained on one experiment predict the other.</b> An additive model learned in Exp1 predicts expression of <b>5,555 genes</b> in Exp2, and the reverse direction predicts 4,006 genes (BH q&lt;0.05). With permuted TF labels the counts are 9 and 0. This is the direct test Reviewer 1 asked for (M6).</li>
<li><b>Combinatorial effects exist but are second-order.</b> A universal test of all 253 TF&times;TF interaction terms finds <b>1,159 reproducible interactions</b> (743 genes; permutation FDR 0.24%). 86% are super-additive. They survive adjustment for each TF's own dose and for load-dependent TF effects (89% retained). Adding them improves cross-experiment prediction only modestly: 5,690 vs 5,555 genes, and 4,229 vs 4,006.</li>
<li><b>NAC4 is the reproducible hub, not ERF5.</b> NAC4 has 1,049 consensus targets (Exp1 vs Exp2 agreement r=0.38), CRF4 has 147 (r=0.24), and ERF5 has 99 (r=&minus;0.07). ERF5 has 1,071 q&lt;0.05 targets in Exp1 but only 30 in Exp2, and its Exp1 signature does not transfer. <b>ERF5&ndash;NAC4 reproduces</b> (59 interactions), but <b>ERF5&ndash;CRF4 does not</b> in the single-cell data (4 interactions, 7 co-targets).</li>
<li><b>CART as currently run is not reproducible at the level of individual trees.</b> With the CART-CELLS settings run independently on each experiment, the root TF (TF1) of a gene's tree is the same in Exp1 and Exp2 for 5.2% of genes, against 5.0% expected by chance. This holds even for genes whose TF1 split has p&lt;1e-10 in both experiments (4.8% vs 4.9%). Trees on covariate-corrected data do predict the other experiment, but less well than the linear model, and their individual split choices remain unstable.</li>
</ul></div>

<h2>1. Data and TF calling</h2>
<p>We kept genes detected in &ge;5% of cells in <i>both</i> experiments (11,208 genes). The endogenous loci of the 27 transfected TFs were excluded as targets, because transgene reads map to them (for example, BEE2 transgene vs locus &rho;=0.78 in Exp2). We kept the 23 GR-TFs detected in &ge;100 cells in both experiments. FBH4, GATA17L, HSFB2A and HYH were dropped because they are essentially absent (&le;94 cells), which matches the manuscript's "23 TFs".
Exp1 cells carry a median of 14 detected GR-TFs and Exp2 cells a median of 10. Cells with no GR-TF read make up 4.9% of Exp1 and 12.7% of Exp2.</p>
<figure>{img('fig1a_tf_detection','48%')} {img('fig1b_tfs_per_cell','50%')}<figcaption>Fig 1. A: fraction of cells with &ge;1 junction read per TF. B: number of TFs detected per cell (Reviewer 1 S5, Reviewer 2 #1/#4).</figcaption></figure>
<p><b>Transfection load is the main confounder.</b> RFP reads and total GR-TF reads per cell are strongly anti-correlated with endogenous UMIs: log RFP vs log gene UMI gives r=&minus;0.68 in Exp1 and &minus;0.63 in Exp2. Load also differs between cell types; for example, Exp1 cortex cells carry a mean of 5.3 TFs and xylem cells 15.0. Every model therefore includes log UMI, log genes detected, log RFP, log total GR-TF reads and library (Exp2). It also includes nine marker-based root cell-type scores and cubic terms of the load variables (see section 3).</p>
<p><b>The matrices look SCTransform-corrected.</b> All values are integers, but per-cell totals vary little (CV 0.18&ndash;0.21). Junction-read counts are therefore not a clean measure of dose, which probably explains why intermediate thresholds behave poorly (Fig 2).</p>

<h2>2. Calibrating the analysis on reproducibility</h2>
<p>For each analysis choice, we ran the same model in Exp1 and Exp2 and asked how well the TF&times;gene effects agree. We also built a null in which each cell's TF profile is permuted among cells of the same library and the same GR-TF-load decile. That permutation keeps the load structure and breaks only the link between TF identity and expression. The choices we kept are those that maximise agreement above this null.</p>
<figure>{img('fig2_calibration_tf_calling','90%')}<figcaption>Fig 2. Correlation of Exp1 and Exp2 TF&times;gene t-statistics for six TF encodings and three covariate sets; crosses show the permutation null (Reviewer 3 #3: sensitivity across thresholds).</figcaption></figure>
<ul>
<li><b>TF calling:</b> presence with &ge;1 junction read, the manuscript's 0.5 threshold, gives the most reproducible binary call (r=0.128 vs null 0.018). At &ge;2 reads agreement falls to the null level (0.020 vs 0.021); at &ge;3 and &ge;5 reads it partly recovers (0.084 and 0.088). The continuous "share of the cell's GR-TF reads" encoding is as good (0.138) and is used as a sensitivity analysis: it gives 2,503 consensus edges at FDR 0.002%.</li>
<li><b>Cell identity:</b> adding marker-based cell-type scores barely changes agreement (0.128 to 0.118 at &ge;1 read), so the TF effects are not explained by cell type (Reviewer 2 #6, Reviewer 3 #5). Data-driven expression PCs, by contrast, remove almost all reproducible signal (5 PCs: r 0.128 to 0.033). The reason is that the PCs carry TF programs: the TF design explains 25% of PC2 in Exp1, driven mainly by ERF5, NAC4 and WRKY18, and 19% of PC3 in Exp2. PCs are therefore not used.</li>
<li><b>Load and interactions:</b> with only linear load terms, the interaction test is badly inflated. Permuted data then give 3,000&ndash;8,000 "consensus" interactions, because TF products act as non-linear load terms. Adding squared and cubic terms of log UMI, log RFP and log GR-TF reads, plus their product, brings the null back to the nominal rate (about 5% of p&lt;0.05) and to about 3 consensus interactions. Polynomials of the number of TFs detected were tested and dropped: they absorb true signal without improving the null.</li>
<li><b>Robust SEs:</b> HC3, the leverage-corrected sandwich estimator. Interactions are only tested when each of the ++, +&minus; and &minus;+ groups has &ge;30 cells in both experiments.</li>
</ul>

<h2>3. Combining the experiments: consensus network</h2>
<p>For each gene we model log1p(CP10k) as a function of the 23 presence indicators plus the covariates, separately in each experiment. Interaction models add one TF<sub>a</sub>&times;TF<sub>b</sub> product at a time, for all 253 pairs. The two experiments are combined with a fixed-effect inverse-variance meta-analysis. An edge is called <b>consensus</b> when it meets three conditions: meta-analysis BH q&lt;0.05, nominal p&lt;0.05 in <i>each</i> experiment, and the same sign in both. Cochran's Q flags heterogeneous effects.</p>
<figure>{img('fig4_permutation_null','90%')}<figcaption>Fig 3. Observed consensus counts compared with 50 permutations of TF labels within library &times; load strata (Reviewer 2 #5).</figcaption></figure>
<div class="tw"><table><tr><th>TF</th><th>TF+ cells Exp1</th><th>TF+ cells Exp2</th><th>q&lt;0.05 Exp1</th><th>q&lt;0.05 Exp2</th><th>r(Exp1,Exp2)</th><th>Consensus targets</th><th>up/down</th></tr>{tf_rows}</table></div>
<figure>{img('fig3a_consensus_targets_per_tf','49%')} {img('fig3b_per_tf_agreement','49%')}</figure>
<figure>{img('fig8_effects_exp1_vs_exp2')}<figcaption>Fig 4. Effect of the TF in Exp1 (x) and Exp2 (y) for genes significant in either experiment. NAC4 and CRF4 fall on the diagonal. TGA1 effects are mostly Exp2-specific and ERF5 effects mostly Exp1-specific.</figcaption></figure>
<div class="box warn"><b>ERF5 needs a careful look before resubmission.</b> In Exp1, ERF5 is the strongest TF:
<ul><li>it has 1,071 q&lt;0.05 targets;</li>
<li>it is TF1 of 2,199 CART-CELLS trees;</li>
<li>it drives expression PC2 (t=&minus;15).</li></ul>
None of this replicates in Exp2, where ERF5 has 30 targets. A signature built from the Exp1 ERF5 targets shows no dose response to ERF5 reads in Exp2 (<code>dose_response_signatures.tsv</code>). The Fig 6 gene AT1G65845 illustrates the problem: its ERF5 effect is &beta;=0.31 (q=1e-4) in Exp1 and &minus;0.05 (ns) in Exp2, with heterogeneity p&lt;1e-4. Its NAC4 effect, by contrast, is consensus. ERF5 also behaves differently in Exp1: it is detected in fewer cells (35% vs 52%), it more often dominates its cell's GR-TF reads (10% vs 2% of ERF5+ cells), and it co-occurs less with other TFs. Did the ERF5 plasmid, prep or amount differ between the two transfections? The DoubleTARGET bulk RNA-seq remains independent support for ERF5&ndash;CRF4, but the single-cell data alone do not reproduce it.</div>

<h2>4. TF&times;TF interactions (universal ANOVA-style test)</h2>
<p>This addresses Reviewer 1 M2 and Reviewer 2 #5. All 253 pairs were tested against all genes, with every main effect and covariate in the model. The terms are used as follows:
<ul>
<li><i>Combinatorial</i> means two TFs regulate the same gene (additive co-targets).</li>
<li><i>Interaction</i> means the effect of the pair differs from the sum of the two single effects on the log scale.</li>
<li><i>Super-additive</i> means the interaction has the same sign as the summed main effects; <i>sub-additive</i> means the opposite sign.</li>
</ul>
Neither type implies physical cooperation.</p>
<p><b>1,159 consensus interactions</b> (null mean 2.8, max 10 in 50 permutations): 1,000 super-additive and 159 sub-additive. Exp1 alone is under-powered for this test. Only 231 interactions reach q&lt;0.05 there, and they replicate at chance level: 4.8% at p&lt;0.05, with 49% sign agreement. The combinatorial calls in the current manuscript, which are based on Exp1 only, should therefore be replaced by these combined calls.</p>
<figure>{img('fig7_interaction_matrix','85%')}<figcaption>Fig 5. Number of reproducible interactions per TF pair. This is the TF-by-TF matrix Reviewer 2 (minor 5) asked for.</figcaption></figure>
<div class="tw"><table><tr><th>Pair</th><th>++ cells Exp1/Exp2</th><th>q&lt;0.05 Exp1</th><th>q&lt;0.05 Exp2</th><th>Consensus interactions</th><th>Consensus co-targets</th></tr>{pair_rows}</table></div>
<p><b>Pairs highlighted in the manuscript</b> (Sup Table 2):</p>
<div class="tw"><table><tr><th>Pair</th><th>++ cells Exp1/Exp2</th><th>Consensus interactions</th><th>Consensus co-targets</th></tr>{paper_rows}</table></div>
<p><b>Robustness.</b> For the 10 pairs with the most interactions, plus ERF5&ndash;CRF4, we refitted with two additional sets of terms. The first is each TF's own dose (log1p reads), which addresses Reviewer 3 #1, the concern that a partner dilutes the dose. The second lets each TF's effect scale with GR-TF load. Of the 438 interactions, 395 were retained with dose terms and 389 with dose plus load modulation. The interaction coefficients correlate with the originals at 0.98&ndash;0.999 after dose adjustment and at 0.85&ndash;0.96 after adding load modulation (<code>interaction_dose_sensitivity.tsv</code>).</p>

<h2>5. Does a model learned in one experiment predict the other?</h2>
<p>This addresses Reviewer 1 M6 and Reviewer 3 #4. Each model is trained on all cells of one experiment and predicts every gene in the other. The score is the correlation between predicted and observed expression after removing the test experiment's covariates. All tuning parameters were chosen on this transfer score:</p>
<ul>
<li>the CART grid covered depth 1&ndash;5, minimum leaf 25/50/100/200 cells or 10% of cells, and cp 0&ndash;0.01;</li>
<li>the ridge penalty for the linear models was chosen the same way.</li>
</ul>
<figure>{img('fig5_cross_experiment_prediction','95%')}<figcaption>Fig 6. Genes significantly predicted in the held-out experiment (BH q&lt;0.05). With permuted test labels at most 35 genes reach significance.</figcaption></figure>
<ul>
<li>The best CART is depth 3 with a 10%-of-cells minimum leaf and no pruning, which is the CART-CELLS architecture, but only when it is fitted on covariate-corrected expression. Depth 2&ndash;5 perform alike and a single split does almost as well. Deeper combinatorial structure adds little that transfers.</li>
<li>The linear additive model predicts 63&ndash;84% more genes than the best tree.</li>
<li>Adding interactions helps only when they are heavily shrunk: 5,690 vs 5,555 genes and 4,229 vs 4,006. Lightly penalised interaction models predict worse (overfitting).</li>
<li>Tree structure is unstable even when prediction works. The root TF of the tuned trees agrees between experiments for 6.4% of genes, against 6.2% by chance. Among the 1,319 genes predicted in both directions it agrees for 26%.</li>
</ul>
<figure>{img('fig5b_cart_cells_tf1_agreement','80%')}<figcaption>Fig 7. CART-CELLS run as published, with classification trees on counts, GR-TF counts as features, depth 3, 10% leaves and Gini, independently on each experiment. TF1 agreement between experiments equals chance at every significance level.</figcaption></figure>

<h2>6. Power: does the network saturate?</h2>
<figure>{img('fig6_power_vs_cells','75%')}<figcaption>Fig 8. Consensus edges when one experiment is down-sampled, 3 random draws per size (Reviewer 2 #4/#8).</figcaption></figure>
<p>Reproducible edges rise roughly linearly with cell number, and neither curve plateaus. This contradicts the manuscript's split-half "saturation" argument. For the design, more cells, or replicate libraries, will keep adding targets.</p>

<h2>7. Reviewer points this answers</h2>
<div class="tw"><table><tr><th>Point</th><th>What the analysis now provides</th></tr>
<tr><td>R2 #4, R3 #2 (biological replicate)</td><td>Two independent experiments analysed separately, combined by meta-analysis, and replication required for every edge (sections 3&ndash;4).</td></tr>
<tr><td>R1 M6 (does CART predict the 2nd set)</td><td>Train-on-one, predict-the-other for CART and linear models (section 5, Fig 6).</td></tr>
<tr><td>R1 M2, R2 #5 (interaction vs combinatorial)</td><td>Explicit TF<sub>a</sub>&times;TF<sub>b</sub> term for all 253 pairs, with a permutation null, dose and load-modulation controls, and clear definitions (section 4).</td></tr>
<tr><td>R1 M3, R2 #1/#5, R3 #5 (TF load, depth)</td><td>Load, RFP, depth and genes detected are in every model, with non-linear terms; a stratified permutation keeps the load structure (sections 1&ndash;3).</td></tr>
<tr><td>R1 M3(2) (CART favours low-prevalence TFs?)</td><td>Partly confirmed. In Exp1, how often a TF is chosen as TF1 rises as its prevalence moves toward the middle of the range (Spearman &rho;=&minus;0.54 with |prevalence&minus;50%|): ERF5 (35% of cells) is TF1 for 2,199 genes. This trend is absent in Exp2 (&rho;=&minus;0.09). Per-TF counts are in <code>cart_faithful_TF1_counts.tsv</code>. The linear model does not choose among TFs, so it avoids this bias.</td></tr>
<tr><td>R2 #6, R3 #5 (cell identity)</td><td>Marker-based cell types: <code>Exp*_TF_by_celltype.tsv</code> gives TF frequency by type. Adjusting for cell identity leaves the TF effects intact (section 2).</td></tr>
<tr><td>R3 #3 (threshold sensitivity)</td><td>Fig 2: six encodings &times; three covariate sets, each against its null.</td></tr>
<tr><td>R3 #4 (hyperparameters, CV)</td><td>Every parameter chosen on cross-experiment transfer (section 5); grids are in <code>cart_grid_transfer.tsv</code> and <code>transfer_linear_lambda_grid.tsv</code>.</td></tr>
<tr><td>R2 #4/#8 (cells per TF and combination, power)</td><td>Cells per TF in the per-TF table; ++/+&minus;/&minus;+/&minus;&minus; cells per pair in <code>per_pair_summary.tsv</code>; down-sampling in Fig 8.</td></tr>
<tr><td>R3 minor 7 (reproducible environment)</td><td>Scripts with fixed seeds; versions listed below.</td></tr>
</table></div>

<h2>8. Draft Methods text</h2>
<p style="background:var(--soft);padding:12px 16px"><i>Statistical analysis of scTARGET-seq.</i> Two independent scTARGET-seq experiments were analysed: Exp1 (2,313 cells) and Exp2 (6,374 cells from six 10x libraries). Each came from a separate protoplast preparation and transfection. Genes detected in &ge;5% of cells in both experiments (n=11,208) were analysed, excluding the endogenous loci of the transfected TFs. Each GR-TF was called present in a cell when &ge;1 fusion-junction read was detected. This threshold maximised the agreement of TF effects between experiments among the thresholds tested (&ge;1, 2, 3, 5 reads, log dose, and share of GR-TF reads). Expression was modelled per gene as log1p(counts per 10,000), using ordinary least squares with HC3 heteroskedasticity-robust standard errors. The predictors were the presence of all 23 GR-TFs. The covariates were log UMI, log genes detected, log RFP reads, log total GR-TF reads (with squared and cubic terms of these load variables and their product), library, and nine marker-based root cell-type scores. TF&times;TF interactions were tested for all 253 pairs by adding one product term at a time (Frisch&ndash;Waugh&ndash;Lovell), when each single- and double-positive group contained &ge;30 cells in both experiments. Estimates from the two experiments were combined by fixed-effect inverse-variance meta-analysis. An effect was retained when its meta-analysis Benjamini&ndash;Hochberg q was &lt;0.05 and it reached p&lt;0.05 with the same sign in each experiment. Empirical FDRs were obtained by repeating the whole procedure after permuting TF profiles among cells within strata of library &times; GR-TF-load decile (50 permutations). Predictive validity was assessed by training models on one experiment and predicting gene expression in the other. CART hyperparameters (rpart, depth 1&ndash;5, minimum leaf, complexity) and ridge penalties were selected on this cross-experiment criterion.</p>

<h2>9. Open questions for gab</h2>
<ol>
<li><b>Raw UMI counts.</b> Are the raw CellRanger UMI counts for the junction features available? The matrices look SCTransform-corrected, which blurs dose. Raw counts would allow proper count models and a cleaner dose analysis.</li>
<li><b>Exp2 libraries.</b> Are the six barcode suffixes in Exp2 six 10x channels from one transfection? The manuscript's earlier replicate had 6,176 cells; this file has 6,374.</li>
<li><b>ERF5 between experiments.</b> Did the ERF5 construct, prep or amount differ between the two experiments? (Section 3.)</li>
<li><b>Paper's GRN parameters.</b> Which CART filter thresholds produced the manuscript's GRN? They would let us compute its exact reproducibility, rather than the default-parameter version shown here.</li>
</ol>

<h2>Files</h2>
<p><code>results/consensus_ge1/</code> contains:
<ul>
<li><code>main_effects_all_significant.tsv.gz</code> and <code>interactions_all_significant.tsv.gz</code>: per-experiment and meta statistics, with a consensus flag;</li>
<li><code>per_TF_summary.tsv</code> and <code>per_pair_summary.tsv</code>: summaries with cell counts per group;</li>
<li><code>consensus_interactions_classified.tsv</code>;</li>
<li><code>permutation_null.tsv</code>.</li>
</ul>
<code>results/cart/</code> holds the CART-CELLS faithful, grid and tuned outputs. <code>results/consensus_share/</code> holds the sensitivity analysis with the continuous encoding. Software: Python 3.13.12 (numpy 2.4.4, scipy 1.18.1, pandas 3.0.6, statsmodels 0.15.0) and R 4.5.2 (rpart 4.1.24, ggplot2 4.0.1, data.table 1.17.8). Scripts 01&ndash;12 run in order; seeds are set in each script.</p>
</body></html>"""
open(os.path.join(BASE, "REPORT_reproducibility.html"), "w").write(html)
print("report written", len(html) // 1024, "KB")
