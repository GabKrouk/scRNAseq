#!/usr/bin/env python3
"""Step 32: QC from the raw CellRanger filtered matrices (Reviewer 2 #9, Reviewer 3 minor 6).
Per cell: UMI, genes, % mitochondrial (ATMG), % plastid (ATCG), GR-TF / RFP / GFP reads
(Exp2 only; the Exp1 CellRanger 4.0.0 run has no transgene features), match to the
cells of the normalised matrices used in the paper, and MAD-based QC flags."""
import os, sys, gzip
import numpy as np, pandas as pd, scipy.io as sio
sys.path.insert(0, os.path.dirname(__file__))
from common import RES, DATA
TOP = "/Users/gabkrouk/Library/CloudStorage/Dropbox/Science/projects/analyse_gab_scTARGET"
RAWD = os.path.join(TOP, "raw_data")
OUT = os.path.join(RES, "qc_raw"); os.makedirs(OUT, exist_ok=True)
summary = {}
for e in ("Exp1", "Exp2"):
    M = sio.mmread(os.path.join(RAWD, f"{e}_matrix.mtx")).tocsr()
    feat = pd.read_csv(os.path.join(TOP, f"{e}_features.tsv.gz"), sep="\t", header=None, names=["id", "name", "type"])
    bc = pd.read_csv(os.path.join(TOP, f"{e}_barcodes.tsv.gz"), header=None)[0]
    isgene = feat.id.str.match(r"^AT[1-5CM]G\d{5}") | (feat.name == "atp4")     # atp4 = mitochondrial gene with a non-AGI id
    mito = (feat.id.str.startswith("ATMG") | (feat.name == "atp4")).values; plast = feat.id.str.startswith("ATCG").values
    nuc = (isgene & ~feat.id.str.startswith(("ATMG", "ATCG"))).values
    G = M[isgene.values]
    tot = np.asarray(G.sum(0)).ravel()
    df = pd.DataFrame({"barcode": bc.values,
                       "umi": tot, "genes": np.asarray((G > 0).sum(0)).ravel(),
                       "pct_mito": 100 * np.asarray(M[mito].sum(0)).ravel() / tot,
                       "pct_plastid": 100 * np.asarray(M[plast].sum(0)).ravel() / tot})
    extra = feat[~isgene]
    if len(extra):
        for _, r in extra.iterrows():
            df[r["name"]] = np.asarray(M[r.name].todense()).ravel()
        tfcols = [c for c in extra.name if c not in ("RFP", "GFP")]
        df["tf_reads"] = df[tfcols].sum(axis=1); df["n_tf"] = (df[tfcols] > 0).sum(axis=1)
    # match to normalised matrix cells (".k" suffix there, "-k" here)
    norm_cells = pd.read_csv(os.path.join(DATA, f"{e}_cells.tsv"), sep="\t").cell.str.replace(".", "-", regex=False)
    df["in_paper_matrix"] = df.barcode.isin(set(norm_cells))
    # MAD-based flags (log scale), computed per library
    df["library"] = df.barcode.str.split("-").str[-1]
    def mad_flag(x, nmad=3, side="low"):
        lx = np.log10(x); med = np.median(lx); mad = 1.4826 * np.median(np.abs(lx - med))
        return lx < med - nmad * mad if side == "low" else lx > med + nmad * mad
    flags = []
    for lib, g in df.groupby("library"):
        f = pd.DataFrame(index=g.index)
        f["low_umi"] = mad_flag(g.umi); f["low_genes"] = mad_flag(g.genes)
        f["high_umi"] = mad_flag(g.umi, side="high")
        pm = g.pct_mito; f["high_mito"] = pm > np.median(pm) + 3 * 1.4826 * np.median(np.abs(pm - np.median(pm)))
        f["high_plastid_5pct"] = g.pct_plastid > 5
        flags.append(f)
    df = df.join(pd.concat(flags))
    df.to_csv(os.path.join(OUT, f"{e}_cell_qc.tsv.gz"), sep="\t", index=False)
    s = dict(cells_cellranger=len(df), cells_in_paper_matrix=int(df.in_paper_matrix.sum()),
             paper_matrix_cells=len(norm_cells), paper_cells_not_in_cellranger=int((~norm_cells.isin(set(df.barcode))).sum()),
             umi_median=float(df.umi.median()), genes_median=float(df.genes.median()),
             pct_mito_median=float(df.pct_mito.median()), pct_mito_p95=float(df.pct_mito.quantile(.95)),
             pct_plastid_median=float(df.pct_plastid.median()), pct_plastid_over5=int((df.pct_plastid > 5).sum()))
    for c in ("low_umi", "low_genes", "high_umi", "high_mito", "high_plastid_5pct"):
        s[f"flag_{c}"] = int(df[c].sum())
        s[f"flag_{c}_in_paper"] = int((df[c] & df.in_paper_matrix).sum())
    anyflag = df[["low_umi", "low_genes", "high_umi", "high_mito", "high_plastid_5pct"]].any(axis=1)
    s["any_flag"] = int(anyflag.sum()); s["any_flag_in_paper"] = int((anyflag & df.in_paper_matrix).sum())
    if "tf_reads" in df:
        s["frac_cells_no_TF_reads"] = float((df.tf_reads == 0).mean()); s["RFP_pos"] = float((df.RFP > 0).mean()); s["GFP_pos"] = float((df.GFP > 0).mean())
    summary[e] = s
    print(e, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in s.items()})
pd.DataFrame(summary).to_csv(os.path.join(OUT, "qc_summary.tsv"), sep="\t")
