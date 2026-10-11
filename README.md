# scRNAseq — analysis code for single-cell TARGET-seq (scTARGET-seq)

Code for *Combinatorial effects of transcription factors revealed by single-cell TARGET-seq*
(Emanuel, Carluer, Hinckley et al.; revision of The Plant Cell TPC-2026-0867).
This repository contains code only; no data.

## Data needed
- CellRanger 7.1 matrices of Exp1 and Exp2 (`Exp{1,2}_matrix.mtx.gz`, `_features.tsv.gz`, `_barcodes.tsv.gz`), INRAE Dataverse https://doi.org/10.57745/SUIM1I
- Bulk ERF5/CRF4 DoubleTARGET counts (`counts.txt`) and metadata
- TAIR GO annotation (`tair.gaf.gz`) and `go-basic.obo.gz` for GO enrichment

Input paths are set at the top of each script (some are absolute and must be edited).

## Layout
- `analysis/scripts/` — numbered steps, run in order.
  - `01` raw counts, cell QC, doublets
  - `02–12` TF detection, calibration of TF encoding and covariates, cell-type scores, consensus models (`05`), CART comparison (`06`), cross-experiment prediction (`07`), power (`09`)
  - `13–31` four-state summaries, qPCR genes, bulk–single-cell bridge, shared dynamics
  - `32–37` QC on raw matrices and SCTransform artefact
  - `38–46` non-additive interactions (load-modulated dose models), robustness, held-out test, Supplemental Data Set 1
  - `47` final model with 50 permutations
  - `48–64` paper figures, bulk FDR re-analysis, GO, cell-type robustness, stress markers, network

## Environment
Python 3.13 (numpy 2.4, scipy 1.18, pandas 3.0, statsmodels 0.15) and
R 4.5.2 (rpart 4.1.24, ggplot2 4.0.1, data.table 1.17.8, gridExtra, jsonlite). Random seeds are fixed in each script.
