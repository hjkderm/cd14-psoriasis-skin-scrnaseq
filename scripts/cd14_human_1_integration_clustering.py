"""
Human skin scRNA-seq (GEO-derived): QC -> HVG selection -> scVI integration -> Leiden clustering

Pipeline overview:
  I.   Load per-sample h5ad files, QC filter, merge
  II.  Highly variable gene (HVG) selection (with manual marker rescue)
  III. Integration with scVI
  IV.  Leiden clustering on the scVI latent space + DEG extraction
  V.   Level-0 cell type annotation
"""

from warnings import filterwarnings
filterwarnings('ignore')

import os
import pandas as pd
import numpy as np

import scanpy as sc

from tqdm import tqdm
from glob import glob

import random
random.seed(0)

# ============================================================
# I. Load samples, QC filter, merge
# ============================================================

# ---- Sample list (from metadata) ----
# meta.csv contains one row per sample; the 'SerSam' column holds the sample
# identifier used to name each per-sample h5ad file.
meta = pd.read_csv('data/meta.csv')
samples = meta['SerSam'].tolist()

raw_h5ad_dir = 'data/01_h5ad_raw/'


alist = []

for sam in tqdm(samples):
    try:
        loc = f'{raw_h5ad_dir}{sam}.h5ad'

        adata = sc.read(loc)
        adata.var_names_make_unique()

        adata.obs.index = [sam + '@' + x for x in adata.obs.index]
        adata.obs['SerSam'] = sam

        # QC metrics
        adata.var['mt'] = adata.var_names.str.startswith("MT-")
        adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
        adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

        sc.pp.calculate_qc_metrics(
            adata, qc_vars=['mt', "ribo", "hb"],
            percent_top=None, log1p=False, inplace=True
        )

        sc.pp.scrublet(adata, verbose=False)

        adata = adata[
            (adata.obs['total_counts'] > 1000)
            & (adata.obs['n_genes_by_counts'] > 500)
            & (adata.obs['n_genes_by_counts'] < 5000)
            & (adata.obs['doublet_score'] < 0.2)
        ].copy()

        adata.layers['count'] = adata.X.copy()

        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.raw = adata.copy()
        adata.write(f'{outdir}00_all/{sam}.h5ad')
        alist.append(adata)

    except Exception as e:
        print(f"Failed to process {sam}: {e}")

merged = sc.concat(alist)

# remove predicted doublets
merged = merged[~merged.obs['predicted_doublet']]

# ============================================================
# II. Highly variable genes (HVGs)
# ============================================================

sc.pp.highly_variable_genes(
    merged,
    subset=False,  # keep full gene set for now; check markers first, subset later
    layer='count',
    flavor='seurat_v3',
    n_top_genes=2000,
    span=0.5,
    batch_key="SerSam",
)

# ---- Check specific genes are included in the HVGs ----
hvg_genes = merged.var['highly_variable']
print(hvg_genes.sum())

marker_genes = [
    # Myeloid
    'LYVE1', 'CD74', 'LYZ', 'CLEC9A', 'GSN', 'MNDA', 'LAMP3', 'FCGR2A', 'C1QA',
    'RAB32', 'CD68', 'LST1', 'CXCL8', 'ITGAX', 'CLEC7A', 'CTSZ', 'TYROBP',
    'SPI1', 'FGL2', 'PTGS2', 'CTSS', 'CSTA',
    'CD207', 'MRC1',
    'KYNU', 'CD163', 'CD14',

    # TNK
    'PTPRC', 'KLRD1', 'CD2', 'CST7', 'XCL1', 'CD3D', 'LTB', 'STAT4', 'CSF2',
    'TNFRSF9', 'GZMB', 'IL7R', 'PI16', 'CD69',
    'CXCR4',
    'CCL5',
    "NKG7", "GNLY",
    'CD7',  

    # B
    'IGHV3-7', 'KMO', 'SPIB', 'CD79A', 'SLAMF7', 'IGHG4',
    'IGKC', 'IGHA1', 'DERL3', 'IGHG1', 'JCHAIN', 'IGHG3', 'IGHM',
    'CD19',
    'MS4A1',
    'IGHV1-18', 'JSRP1', 'IGLV3-1',

    # Keratinocytes
    'SYT8', 'SERPINB5', 'PERP', 'DSC1', 'DST', 'CALML3',
    'SPRR2D', 'GJB6', 'CHI3L1', 'DSG1', 'LGALS7B', 'SBSN', 'KRT1', 'PKP1',
    'CALML5', 'S100A14', 'KRTDAP', 'LCE1C', 'CRABP2', 'DSP', 'NRARP',
    'SAA1', 'LAMB3', 'KRT6B', 'KRT6C', 'DMKN', 'COL17A1', 'KRT10', 'KRT16', 'KRT17',
    'KRT14', 'DEFB1', 'KRT23', 'DSC3',
    'KRT5', 'KRT2', 'EPCAM',
    'KRT9',

    # EccrineG_sec
    'S100A1',
    'KRT19', 'CSRP1', 'CLDN3', 'SNORC', 'KRT7', 'PPP1R1B', 'CA6', 'AZGP1',
    'SCGB1D2', 'KRT18', 'SCGB1B2P', 'SCGB2A2', 'KRT15',
    'PLA2G2A', 'MUCL1', 'AKR1C1',
    'DCD',

    # Fibroblasts
    'APOD',
    'LUM', 'DCN', 'COL1A1', 'COL6A2', 'RBP4', 'PLAC9', 'COCH',
    'CCDC80', 'COL1A2', 'C7', 'SCARA5',
    'COL3A1',
    'APCDD1',  
    'ASPN', 'POSTN', 'TNN', 'COL11A1',  
    'RARRES2',

    # SmoothM
    'PLN', 'C11orf96', 'MUSTN1', 'TAGLN', 'MCAM', 'DES', 'STEAP4', 'ACTG2',
    'EGFL6', 'LMOD1', 'AEBP1', 'MYL9',
    "MYLK", "ACTA2",

    # Endothelial_Vascular
    'PTPRB', 'GNG11', 'PCAT19', 'IFI27', 'AQP1', 'SPARCL1', 'ADGRL4',
    'DNASE1L3', 'CD93', 'RBP7', 'EMCN', 'CX3CL1', 'SOX17', 'ADIRF',
    'PECAM1', 'CD34',
    'VWF',

    # Endothelial_Lymphatic
    'MMRN1', 'PKHD1L1', 'ECSCR', 'TFPI', 'AKAP12', 'RAMP2', 'ANGPT2',
    'CAVIN2', 'EFEMP1',
    'LYVE1',

    # MastCell
    'CPA3', 'TPSAB1', 'TPSB2', 'HPGDS', 'CTSG', 'KIT',
    'GATA2',

    # Melanocyte
    'PMEL', 'DCT', 'TYRP1', 'MLANA', 'TRPM1',
    'CAPN3',

    # Nerve
    'PLP1', 'NRXN1', 'GPM6B',
    'SORBS2', 'SORCS1',
    "PMP22",
]

omitted_markers = [gene for gene in marker_genes
                    if gene not in merged.var_names[hvg_genes].tolist()]
print(f"Omitted markers in HVGs: {omitted_markers}")

# ---- Manually add omitted markers to HVGs ----
for gene in omitted_markers:
    if gene in merged.var_names:
        merged.var.loc[gene, 'highly_variable'] = True

updated_hvg_genes = merged.var['highly_variable']
print(f"Total number of HVGs after appending: {updated_hvg_genes.sum()}")

# ---- Subset to HVGs ----
merged_hvg = merged[:, merged.var['highly_variable']].copy()

# ============================================================
# III. Integration: scVI
# ============================================================

import torch
import scvi

scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)
torch.set_float32_matmul_precision("high")

# ---- Set up scVI model with HVGs ----
scvi.model.SCVI.setup_anndata(
    merged_hvg,
    layer='count',
    batch_key='SerSam',
    continuous_covariate_keys=["pct_counts_mt", "pct_counts_ribo"],
)

# ---- Create and train the model ----
model = scvi.model.SCVI(merged_hvg)
model.train()

model_dir = 'models/scvi_model_human_skin/'
model.save(model_dir)


# ---- Obtain model outputs ----
SCVI_LATENT_KEY = "X_scVI"
latent = model.get_latent_representation()
merged_hvg.obsm[SCVI_LATENT_KEY] = latent

SCVI_NORMALIZED_KEY = "scvi_normalized"
merged_hvg.layers[SCVI_NORMALIZED_KEY] = model.get_normalized_expression(library_size=10e4)


import matplotlib.pyplot as plt
import seaborn as sns

# ---- Visualization with batch correction (scVI latent space) ----
sc.pp.neighbors(merged_hvg, use_rep=SCVI_LATENT_KEY)
sc.tl.umap(merged_hvg, min_dist=0.3)

sc.pl.umap(merged_hvg, color=["SerSam"], frameon=False, palette=palette, legend_loc='none')

merged_hvg.write('data/cd14_merged_hvg_human_skin.h5ad')

# ============================================================
# IV. Clustering on the scVI latent space
# ============================================================
SCVI_CLUSTERS_KEY = "leiden_scVI"

for res in [0.02, 0.1]:
    sc.tl.leiden(merged_hvg, key_added=f"{SCVI_CLUSTERS_KEY}_{res:3.1f}", resolution=res)

sc.pl.umap(merged_hvg, color=["leiden_scVI_0.0", "leiden_scVI_0.1"], legend_loc="on data")

for res in [0.2, 0.5]:
    sc.tl.leiden(merged_hvg, key_added=f"{SCVI_CLUSTERS_KEY}_{res:3.1f}", resolution=res)

sc.pl.umap(merged_hvg, color=['leiden_scVI_0.2', "leiden_scVI_0.5"], legend_loc="on data")


def extract_degs(adata, groupby, output_path, key_added=None,
                  logfc_thresh=1.5, pval_thresh=0.05, numeric_sort=True):
    """Run rank_genes_groups (Wilcoxon) and save filtered DEGs per cluster to CSV."""
    rank_key = key_added or 'rank_genes_groups'
    sc.tl.rank_genes_groups(adata, groupby=groupby, key_added=rank_key, method="wilcoxon")
    rank_genes_groups = adata.uns[rank_key]
    groups = rank_genes_groups['names'].dtype.names
    sorted_groups = sorted(groups, key=lambda x: int(x)) if numeric_sort else groups

    result_dfs = []
    for group in sorted_groups:
        df = pd.DataFrame({
            'Cluster': group,
            'Gene': rank_genes_groups['names'][group],
            'Logfoldchange': rank_genes_groups['logfoldchanges'][group],
            'P-value': rank_genes_groups['pvals'][group],
            'Adjusted P-value': rank_genes_groups['pvals_adj'][group],
            'Score': rank_genes_groups['scores'][group],
        })
        filtered_df = df[(df['Logfoldchange'] > logfc_thresh) & (df['P-value'] < pval_thresh)]
        result_dfs.append(filtered_df)

    all_degs = pd.concat(result_dfs, ignore_index=True)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    all_degs.to_csv(output_path, index=False)
    return all_degs


os.makedirs('results/cluster_markers_DEG', exist_ok=True)
extract_degs(merged_hvg, "leiden_scVI_0.2", "results/cluster_markers_DEG/cluster_markers_DEG_0.2.csv")

for res in [1.0]:
    sc.tl.leiden(merged_hvg, key_added=f"{SCVI_CLUSTERS_KEY}_{res:3.1f}", resolution=res)

sc.pl.umap(merged_hvg, color=['leiden_scVI_1.0'], legend_loc="on data")

merged_hvg.write('data/merged_hvg_human_skin.h5ad')

# ---- Check outlier and small clusters ----
extract_degs(
    merged_hvg, "leiden_scVI_0.5",
    "results/cluster_markers_DEG/cluster_markers_DEG_0.5.csv",
    key_added='rank_genes_scVI_0.5',
)

merged_hvg.write('data/cd14_merged_hvg_human_skin.h5ad')


# ============================================================
# V. Level-0 cell type annotation
# ============================================================

merged_hvg.obs["cell_type_lv0"] = merged_hvg.obs["leiden_scVI_0.2"].map({
    '0': 'Keratinocytes',
    '3': 'Keratinocytes',
    '14': 'Keratinocytes',

    '1': 'Fibroblasts',
    '2': 'TNK',
    '4': 'Myeloid',
    '5': 'Endothelial_Vascular',
    '6': 'SmoothMuscle',
    '7': 'Melanocytes',
    '8': 'B',
    '9': 'Endothelial_Lymphatic',
    '10': 'EccrineG_sec',
    '11': 'MastCell',
    '12': 'Nerve',

})


sc.pl.umap(merged_hvg, color=["cell_type_lv0"], legend_loc="on data", frameon=False)

# ---- Dotplot & DEG for cluster marker genes ----
marker_genes_lv0 = {
    "Keratinocytes": ["KRT1", "KRT5"],
    "EccrineG_sec": ['KRT7', 'AQP5'],
    "TNK": ["CD3E", "NKG7"],
    "B": ['CD79A', 'JCHAIN'],
    "Myeloid": ['KYNU', 'CD68', 'LYZ'],
    "MastCell": ["TPSB2", "TPSAB1"],
    "Fibroblasts": ["DCN", "CFD"],
    "SmoothM": ["TAGLN", 'ACTA2'],
    "Endothelial_Vascular": ["PECAM1", "ACKR1"],
    "Endothelial_Lymphatic": ['MMRN1', "LYVE1"],
    "Melanocyte": ["PMEL", "MLANA"],
    "Nerve": ["CDH19", 'NRXN1'],
}

sc.pl.dotplot(merged_hvg, marker_genes_lv0, groupby="cell_type_lv0", standard_scale="var")

# ---- Copy lv0 annotation and UMAP embedding back onto the un-subsetted object ----
merged = merged[merged_hvg.obs_names].copy()

assert all(merged.obs_names == merged_hvg.obs_names), "Error: cell indices do not match!"

merged.obs['annot'] = merged_hvg.obs['annot'].copy()

for key in ['X_pca', 'X_umap', 'X_scVI']:
    merged.obsm[key] = merged_hvg.obsm[key].copy()

merged.write('data/cd14_merged_human_skin.h5ad')
