# Single-cell RNA-seq analysis code

Analysis code for the single-cell RNA-seq processing, integration, clustering,
cell-type annotation and cell–cell communication analyses reported in the paper.
It covers two datasets: a **mouse** skin dataset generated in this study, and a
**human** skin dataset assembled from publicly available samples.

Processing, integration, clustering and annotation are done in **Python** (scanpy /
scvi-tools). The annotated mouse object is then handed over to **R** for the
cell–cell communication analysis (CellChat).

---

## Repository contents

```
.
├── README.md
├── requirements.txt
├── LICENSE
├── .gitignore
├── mouse_clustering_annotation.py
├── mouse_cellchat.R
├── mouse_cellchat_per_sample_permutation.R
├── human_integration_clustering.py
├── human_myeloid_subclustering.py
├── human_integration_batch_evaluation.py
└── human_myeloid_composition_analysis.py
```

| Order | Script | Language | Dataset | What it does |
| --- | --- | --- | --- | --- |
| 1 | `mouse_clustering_annotation.py` | Python | Mouse skin | 10x matrices → QC → scVI integration → Leiden clustering → marker-based annotation → per-group cell counts |
| 2 | `mouse_cellchat.R` | R | Mouse skin | Cell–cell communication analysis comparing *Cd14*-KO and WT (CellChat) |
| 3 | `mouse_cellchat_per_sample_permutation.R` | R | Mouse skin | Per-mouse CellChat with an exact permutation test and direction-concordance metric |
| 4 | `human_integration_clustering.py` | Python | Human skin | Per-sample QC → HVG selection → scVI integration → Leiden clustering → level-0 (lv0) cell-type annotation |
| 5 | `human_myeloid_subclustering.py` | Python | Human skin, myeloid | Subset to myeloid cells → Harmony re-integration → Leiden sub-clustering → myeloid subtype annotation (lv1) |
| 6 | `human_integration_batch_evaluation.py` | Python | Human skin | Batch-effect evaluation: before/after UMAPs and iLISI/cLISI metrics for both integrations |
| 7 | `human_myeloid_composition_analysis.py` | Python | Human skin, myeloid | Sample-size downsampling, per-dataset composition breakdown, CFLAR expression proportion |

The mouse and human analyses are independent of each other. Within each, the
scripts must be run in the order shown, since the later scripts read the objects
the earlier ones write.

> **Python → R hand-off.** The R scripts do not read the `.h5ad` output directly.
> Convert the annotated mouse object to a Seurat object saved as
> `data/mouse_skin_for_cellchat.rds` first (for example with `sceasy`,
> `zellkonverter`, or by exporting the counts matrix and the `annot` / `group`
> columns of `.obs`). The R scripts expect the cell-type labels in `annot`, the
> genotype in `group` (`KO` / `WT`), and the mouse identifier in `sample`.

These were originally developed as Jupyter notebooks and have been converted to
plain `.py` scripts for readability. They are written to be read top to bottom as
a record of the analysis; several steps (clustering resolution choice, marker
inspection, cluster annotation) were performed interactively, so running the files
end to end unattended is not expected to reproduce every figure without manual
review at the annotation steps.

---

## Data

**Raw and processed data files are not included in this repository.** The scripts
expect data in a local `data/` directory, which is excluded via `.gitignore`.

### Expected input layout

```
data/
├── meta.csv                     # human: sample table; must contain a 'SerSam' column
├── 00_matrix/                   # mouse: 10x CellRanger output, one folder per sample
├── 01_h5ad_raw/                 # per-sample .h5ad files, named <sample>.h5ad
├── 02_h5ad_filtered/            # written by the QC steps
├── mouse_skin_merged_hvg.h5ad   # mouse: HVG subset carrying the scVI latent space
├── mouse_skin_annotated.h5ad    # mouse: full gene set with the cell-type annotation
├── mouse_skin_for_cellchat.rds  # mouse: annotated object converted to Seurat, for CellChat
├── human_skin_merged_hvg.h5ad   # human: HVG subset carrying the scVI latent space
├── human_skin_annotated.h5ad    # human: full gene set with the lv0 annotation
├── human_skin_myeloid.h5ad      # human: myeloid subset with the lv1 annotation
└── human_skin_myeloid_final.h5ad  # human: myeloid subset with the final lv1 annotation
```

Paths in the scripts are relative placeholders (`data/`, `models/`, `results/`).
Adjust them at the top of each script to match your own layout.

---

## Installation

### Python

```bash
# Python 3.9–3.11 recommended
python -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

A GPU is strongly recommended for the scVI steps (the mouse and human integration
scripts); they will run on CPU but slowly. The Harmony step in the myeloid script
runs fine on CPU.

### R

The two `.R` scripts require R:

```r
install.packages(c("Seurat", "patchwork", "ggplot2", "future"))
# CellChat is installed from GitHub
if (!require("devtools")) install.packages("devtools")
devtools::install_github("jinworks/CellChat")
```

### Key package versions used

| Package | Version |
| --- | --- |
| scanpy | 1.10.2 |
| scvi-tools | 1.1.6.post2 |
| scib | 1.1.7 |
| CellChat | 2.1.2 |
| Seurat | 5.3.0 |

See `requirements.txt` for the full Python dependency list.

---

## Analysis overview

### Mouse skin — `mouse_clustering_annotation.py`

Six samples (three *Cd14*-KO, three WT) across three batches.

- **QC**: per-sample filtering that keeps cells with total counts >1,000, 250–6,000
  genes, mitochondrial fraction <5% and a Scrublet doublet score <0.2 (Scrublet as
  implemented in Scanpy). Samples with fewer than 1,000 cells after filtering are
  excluded. A stricter doublet-score cutoff (≤0.1) is applied after merging.
- **Integration**: scVI, with `sample` as batch key and `pct_counts_mt` as a continuous
  covariate. 2,000 HVGs (`seurat_v3`), with a curated marker panel added back manually
  so that lineage markers are not lost during HVG subsetting.
- **Clustering**: Leiden on the scVI latent space at resolutions 0.1 / 0.2 / 0.5 / 1.0,
  with the myeloid compartment re-clustered at higher resolution.
- **Annotation**: manual, based on cluster-level differentially expressed genes
  (Wilcoxon) and a curated marker panel.

### Mouse cell–cell communication — `mouse_cellchat.R`

CellChat analysis of the annotated mouse object, run separately for *Cd14*-KO and WT
and then compared. Uses the mouse CellChat ligand–receptor database. Produces the
differential interaction heatmaps, the information-flow ranking, and the CXCL/CCL
pathway network plots, and writes the full ligand–receptor interaction tables for
each genotype to `results/`.

### Mouse per-sample validation — `mouse_cellchat_per_sample_permutation.R`

Pooling three mice per genotype into one CellChat object collapses the biological
replicates, so genotype differences cannot be tested for reproducibility across
mice. This script instead runs CellChat independently on each of the six mice,
extracts one pathway-strength value per mouse per pathway, and compares genotypes
with an exact permutation test over all C(6,3) = 20 relabelings, with BH
correction. Because n = 3 per group places a hard floor of 0.1 on the two-sided
p-value, a direction-concordance metric across the 3 × 3 pairwise comparisons is
reported alongside.

### Human skin — `human_integration_clustering.py`

Publicly available human skin samples (378,424 cells after QC), listed in
`data/meta.csv` via the `SerSam` column.

- **QC**: per-sample filtering that keeps cells with total counts >1,000, 500–5,000
  genes and a Scrublet doublet score <0.2 (Scrublet as implemented in Scanpy); cells
  classified as predicted doublets are removed after merging.
- **Integration**: scVI, batch key `SerSam`, with `pct_counts_mt` and `pct_counts_ribo`
  as continuous covariates.
- **Clustering**: Leiden at resolutions 0.02 / 0.1 / 0.2 / 0.5 / 1.0.
- **Annotation (lv0)**: 12 broad compartments — keratinocytes, fibroblasts, T/NK,
  myeloid, vascular and lymphatic endothelium, smooth muscle, melanocytes, B cells,
  eccrine gland, mast cells and nerve. Low-quality clusters labelled `OUT` are removed.

### Human myeloid subset — `human_myeloid_subclustering.py`

- **Re-integration**: the myeloid compartment is subset from the lv0 object and
  re-integrated independently with **Harmony** (`SerSam` as batch key), rather than
  scVI, to resolve finer myeloid structure.
- **Clustering**: Leiden at resolutions 0.2 / 0.5 / 1.0 / 2.0 / 3.0; the annotation is
  hierarchical, with the coarse labels taken from resolution 0.2 and individual
  populations refined at higher resolutions.
- **Annotation (lv1)**: cDC1, cDC2, mregDC, moDC, pDC,
  Langerhans cells, M1-like and M2-like macrophages, and neutrophils.

### Human batch-effect evaluation — `human_integration_batch_evaluation.py`

Quantifies how well integration removed study-level batch effects, at both the
whole-tissue (scVI) and myeloid (Harmony) levels. For each, a UMAP computed from
the uncorrected PCA is compared against the integrated embedding, coloured by
dataset, and iLISI (batch mixing) and cLISI (cell-type separation) are computed
with the `scib` package using dataset as the batch variable. For the myeloid
UMAPs, cells per dataset are capped for visualization only, so that the largest
datasets do not visually dominate; all metrics use the full data.

### Human myeloid composition — `human_myeloid_composition_analysis.py`

- **Sample-size downsampling**: the public cohort is imbalanced (83 healthy vs 18
  psoriasis samples). Healthy samples are repeatedly downsampled to match the
  psoriasis sample number over 50 iterations, and myeloid subtype proportions are
  recalculated and compared against the full-cohort values.
- **Per-dataset breakdown**: myeloid subtype composition per dataset, to check
  whether the overall composition is driven by any single dataset.
- **CFLAR expression proportion**: fraction of CFLAR-expressing cells per myeloid
  subtype, reported alongside mean expression because mean expression alone is
  vulnerable to dropout for sparsely detected genes.

---

## License

Released under the MIT License — see [LICENSE](LICENSE).
