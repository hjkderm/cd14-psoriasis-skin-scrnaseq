# Single-cell RNA-seq analysis code

This repository contains the scripts used for single-cell RNA-seq processing,
integration, clustering and cell-type annotation reported in the paper. It covers
two datasets: a **mouse** skin dataset generated in this study, and a **human**
skin dataset assembled from publicly available samples.

---

## Repository contents

```
.
├── README.md
├── requirements.txt
├── LICENSE
├── .gitignore
├── mouse_clustering_annotation.py
├── human_integration_clustering.py
└── human_myeloid_subclustering.py
```

| Order | Script | Dataset | What it does |
| --- | --- | --- | --- |
| 1 | `mouse_clustering_annotation.py` | Mouse skin (this study) | 10x matrices → QC → scVI integration → Leiden clustering → marker-based annotation → per-group cell counts |
| 2 | `human_integration_clustering.py` | Human skin (public) | Per-sample QC → HVG selection → scVI integration → Leiden clustering → level-0 (lv0) cell-type annotation |
| 3 | `human_myeloid_subclustering.py` | Human skin, myeloid subset | Subset to myeloid cells → Harmony re-integration → Leiden sub-clustering → myeloid subtype annotation (lv1), merged with finer lv2 labels |

The two human scripts must be run in order: `human_myeloid_subclustering.py` reads
the annotated object produced by `human_integration_clustering.py`. The mouse script
is independent of the other two.

These were originally developed as Jupyter notebooks and have been converted to
plain `.py` scripts for readability. They are written to be read top to bottom as
a record of the analysis; several steps (clustering resolution choice, marker
inspection, cluster annotation) were performed interactively, so running the files
end to end unattended is not expected to reproduce every figure without manual
review at the annotation steps.

---

## Data availability

Sequencing data are deposited in GEO under accession **\<FILL IN: GSExxxxxx\>**
(see the Supplementary Information of the paper for the full accession list,
including the accessions of the publicly available human datasets re-analysed here).

**Raw and processed data files are not included in this repository.** The scripts
expect data in a local `data/` directory, which is excluded via `.gitignore`.

### Expected input layout

```
data/
├── meta.csv                  # human: sample table; must contain a 'SerSam' column
├── meta_cd14.csv             # (not required — mouse groups are parsed from sample names)
├── 00_matrix/                # mouse: 10x CellRanger output, one folder per sample
├── 01_h5ad_raw/              # per-sample .h5ad files, named <sample>.h5ad
├── 02_h5ad_filtered/         # written by the QC steps
└── subsets_myeloid/          # lv2 myeloid sub-cluster objects (cdc2f/m1/m2.h5ad)
```

Paths in the scripts are relative placeholders (`data/`, `models/`, `results/`).
Adjust them at the top of each script to match your own layout.

---

## Installation

```bash
# Python 3.9–3.11 recommended
python -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

A GPU is strongly recommended for the scVI steps (the mouse and human integration
scripts); they will run on CPU but slowly. The Harmony step in the myeloid script
runs fine on CPU.

### Key package versions used

The analysis was last run on 2024-12-26 with:

| Package | Version |
| --- | --- |
| scanpy | 1.10.2 |
| scvi-tools | 1.1.6.post2 |

See `requirements.txt` for the full dependency list.

---

## Analysis overview

### 1. Mouse skin — `mouse_clustering_annotation.py`

Six samples (three *Cd14*-KO, three WT) across three batches.

- **QC**: per-sample doublet detection (Scrublet), then filtering on total counts
  (>1,000), genes per cell (250–6,000) and mitochondrial fraction (<5%). Samples with
  fewer than 1,000 cells after filtering are excluded. A stricter doublet-score cutoff
  (≤0.1) is applied after merging.
- **Integration**: scVI, with `sample` as batch key and `pct_counts_mt` as a continuous
  covariate. 2,000 HVGs (`seurat_v3`), with a curated marker panel added back manually
  so that lineage markers are not lost during HVG subsetting.
- **Clustering**: Leiden on the scVI latent space at resolutions 0.1 / 0.2 / 0.5 / 1.0,
  with the myeloid compartment re-clustered at higher resolution.
- **Annotation**: manual, based on cluster-level differentially expressed genes
  (Wilcoxon) and a curated marker panel.

### 2. Human skin — `human_integration_clustering.py`

Publicly available human skin samples (~930,000 cells after QC), listed in
`data/meta.csv` via the `SerSam` column.

- **QC**: per-sample Scrublet, filtering on total counts (>1,000), genes per cell
  (500–5,000) and doublet score (<0.2); predicted doublets removed after merging.
- **Integration**: scVI, batch key `SerSam`, with `pct_counts_mt` and `pct_counts_ribo`
  as continuous covariates.
- **Clustering**: Leiden at resolutions 0.02 / 0.1 / 0.2 / 0.5 / 1.0.
- **Annotation (lv0)**: 13 broad compartments — keratinocytes, fibroblasts, T/NK,
  myeloid, vascular and lymphatic endothelium, smooth muscle, melanocytes, B cells,
  eccrine gland, mast cells and nerve. Low-quality clusters are labelled `OUT`.

### 3. Human myeloid subset — `human_myeloid_subclustering.py`

- **Re-integration**: the myeloid compartment is subset from the lv0 object and
  re-integrated independently with **Harmony** (`SerSam` as batch key), rather than
  scVI, to resolve finer myeloid structure.
- **Clustering**: Leiden at resolutions 0.2 / 0.5 / 1.0 / 2.0 / 3.0; the annotation is
  hierarchical, with the coarse labels taken from resolution 0.2 and individual
  populations refined at higher resolutions.
- **Annotation (lv1)**: cDC1, cDC2, mregDC, mregDC(FSCN1-negative), moDC, pDC,
  Langerhans cells, M1-like and M2-like macrophages, and neutrophils.
- **lv2 labels** produced by separate sub-clustering of the cDC2, M1 and M2
  populations are merged back onto the myeloid object at the end.

---

## Citation

If you use this code, please cite the paper:

```bibtex
@article{<FILL IN: citation key>,
  title   = {<FILL IN: paper title>},
  author  = {<FILL IN: author list>},
  journal = {<FILL IN: journal>},
  year    = {<FILL IN: year>},
  doi     = {<FILL IN: DOI>}
}
```

## License

Released under the MIT License — see [LICENSE](LICENSE).

## Contact

Questions about the code or analysis: \<FILL IN: corresponding author name and email\>
