"""
CD14 dataset: QC -> SCVI integration -> Leiden clustering -> marker-based annotation

Pipeline overview:
  I.   Load datasets (10x matrix -> h5ad)
  II.  QC and filtering (doublet score, gene/count thresholds)
  III. Integration with scVI + Leiden clustering at multiple resolutions
  IV.  Marker-based cluster annotation
  V.   Cell count summary per group / sample

Note: Data paths below are relative placeholders (data/, models/, results/).
Adjust them to match your own directory layout before running.
"""

from warnings import filterwarnings
filterwarnings('ignore')

import os
import pandas as pd
import numpy as np

import scanpy as sc
import scanpy.external as sce

from tqdm import tqdm
from glob import glob

# ============================================================
# I. Load datasets
# ============================================================

# ---- 10x matrix -> h5ad ----
data_folder = 'data/00_matrix/'
h5ad_folder = 'data/01_h5ad_raw/'

# Sample name -> source folder mapping
dataset_indexes = {
    'KO_S': 'Sample_KO',
    'KO_0927': '0927_KO',
    'KO_0928': '0928_KO',
    'WT_S': 'Sample_WT',
    'WT_0927': '0927_WT',
    'WT_0928': '0928_WT',
}

for index, dataset_file in dataset_indexes.items():
    adata = sc.read_10x_mtx(f'{data_folder}{dataset_file}')
    adata.var_names_make_unique()
    adata.write_h5ad(f'{h5ad_folder}{index}.h5ad')

# ---- Sample list ----
# Sample names follow the pattern '<GROUP>_<batch>' (e.g. 'KO_0927', 'WT_S'),
# so the group (KO/WT) can be derived directly from the name instead of a
# separate metadata file.
sample_names = list(dataset_indexes.keys())


def sample_to_group(sample_name):
    return sample_name.split('_')[0]

# ============================================================
# II. QC
# ============================================================

# ---- Violin plots (per-sample QC metrics) ----
from matplotlib import pyplot as plt

for sam in tqdm(sample_names):
    loc = f'{h5ad_folder}{sam}.h5ad'
    adata = sc.read(loc)

    adata.var['mt'] = adata.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars='mt', percent_top=None, log1p=False, inplace=True
    )

    plt.figure()
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=True,
    )
    plt.tight_layout()
    plt.close()

# ---- Filter ----
outdir = 'data/02_h5ad_filtered/'
os.makedirs(outdir + '00_STANDARD', exist_ok=True)
os.makedirs(outdir + '01_HIGHQ', exist_ok=True)

# QC thresholds
min_n_counts = 1000
min_n_genes = 250
max_n_genes = 6000
max_p_mito = 5

sc.settings.verbosity = 0
N = 1000  # minimum cell count required to keep a sample in the HIGHQ set

res_out = []
res_in = []
alist = []  # note: repopulated from disk in the next block, so values appended
            # here are not actually used downstream -- kept only so the loop
            # runs without a NameError.

for sam in tqdm(sample_names):
    loc = f'{h5ad_folder}{sam}.h5ad'
    adata = sc.read(loc)

    adata.obs.index = [sam + '@' + x for x in adata.obs.index]
    adata.var['mt'] = adata.var_names.str.startswith("mt-")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars='mt', percent_top=None, log1p=False, inplace=True
    )

    sc.external.pp.scrublet(adata, verbose=False)

    adata = adata[
        (adata.obs['total_counts'] > min_n_counts)
        & (adata.obs['n_genes_by_counts'] > min_n_genes)
        & (adata.obs['n_genes_by_counts'] < max_n_genes)
        & (adata.obs['pct_counts_mt'] < max_p_mito)
    ].copy()

    adata.layers['count'] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()
    sc.write(f'{outdir}/00_STANDARD/{sam}.h5ad', adata)

    bdata = adata[(adata.obs['doublet_score'] < 0.2)].copy()

    if bdata.shape[0] < N:
        res_out.append((sam, bdata.shape[0]))
    else:
        sc.write(f'{outdir}/01_HIGHQ/{sam}.h5ad', bdata)
        res_in.append((sam, bdata.shape[0]))
        alist.append(bdata)

# ---- Merge high-quality samples ----
alist = []

for sam in tqdm(sample_names):
    loc = f'{outdir}/01_HIGHQ/{sam}.h5ad'
    adata = sc.read(loc)

    adata.obs['sample'] = sam
    adata.obs['group'] = sample_to_group(sam)

    alist.append(adata)

merged = sc.concat(alist)

# Additional doublet-score filtering (stricter cutoff than the per-sample step above)
print((merged.obs['doublet_score'] > 0.1).sum())
merged = merged[merged.obs['doublet_score'] <= 0.1].copy()

# ============================================================
# III. Integration: scVI
# ============================================================
# ref. https://docs.scvi-tools.org/en/stable/tutorials/notebooks/quick_start/api_overview.html
# ref. https://docs.scvi-tools.org/en/stable/tutorials/notebooks/quick_start/data_loading.html

# ---- Highly variable genes (HVGs) ----
# Plan:
# 1. calculate HVGs, append additional omitted specific markers
# 2. subset --> scVI model training
# 3. transfer the identified cell types to the un-subsetted object
# 4. train a second model with all genes for DE analysis only
#    (compare the result with that using the subsetted object)

sc.pp.highly_variable_genes(
    merged,
    subset=False,  # keep full gene set for now; check markers first, subset later
    layer='count',
    flavor='seurat_v3',
    n_top_genes=2000,
    span=0.5,
    batch_key="sample",
)

# ---- Check that specific markers of interest are included in the HVGs ----
marker_genes = [
    'Col3a1', 'Pecam1', 'Lyve1',
    'Tinagl1', 'Tagln', 'Mylk',
    'H2-Aa', 'H2-Ab1', 'Itgam', 'Itgax', 'Adgre1', 'Itgae', 'Cd207',
    'Csf1r', 'Cd68', 'Mertk', 'Marco', 'Ly6c1', 'Cd163',
    'Trbc1', 'Trac', 'Cd8a', 'Cd3d',
    'Nkg7', 'Gzmb', 'Gzmk',  # NK, cytotoxic T
    'Tcrd', 'Tcrg',  # gamma-delta T cell
    'Rorc', 'Il17a', 'Il17f', 'Ccr6', 'Il23r', 'Stat3',  # Th17
    'Bst2', 'Siglech', 'Tcf4',  # pDC
    'Krt5', 'Krt1', 'Col17a1',
    'Fcer1a', 'Cd55', 'Gata2', 'Mrgprb1',  # mast cell
    'Mlana', 'Tyrp1',  # melanocyte
]

hvg_genes = merged.var['highly_variable']
selected_markers = [g for g in marker_genes if g in merged.var_names[hvg_genes].tolist()]
print(f"Selected markers included in HVGs: {selected_markers}")

# ---- Append omitted markers so they aren't lost during subsetting ----
omitted_markers = [g for g in marker_genes if g not in merged.var_names[hvg_genes].tolist()]
print(f"Omitted markers in HVGs: {omitted_markers}")
print(len(omitted_markers))

for gene in omitted_markers:
    if gene in merged.var_names:
        merged.var.loc[gene, 'highly_variable'] = True

updated_hvg_genes = merged.var['highly_variable']
print(f"Total number of HVGs after appending: {updated_hvg_genes.sum()}")

still_missed = [g for g in marker_genes if g not in merged.var_names[updated_hvg_genes].tolist()]
print(f"Still missed markers in HVGs: {still_missed}")

merged.write('data/mouse_skin_annotated.h5ad')

# ---- Subset to HVGs ----
merged_hvg = merged[:, merged.var['highly_variable']].copy()
# .copy() needed: the subset must be a fully independent object, not a view
merged_hvg.write('data/mouse_skin_merged_hvg.h5ad')

# ---- Set up scVI model with HVGs ----
merged_hvg = sc.read('data/mouse_skin_merged_hvg.h5ad')

import torch
import seaborn as sns
import scvi

scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)
torch.set_float32_matmul_precision("high")

scvi.model.SCVI.setup_anndata(
    merged_hvg,
    layer='count',
    batch_key='sample',
    continuous_covariate_keys=["pct_counts_mt"],
    # Even with a single covariate, continuous_covariate_keys must be a list
)

# ---- Train the model ----
# ref. https://docs.scvi-tools.org/en/stable/api/reference/scvi.model.SCVI.html
model = scvi.model.SCVI(merged_hvg)
model.train()

model_dir = 'models/scvi_model_cd14/'
model.save(model_dir)

# To reload later:
# model = scvi.model.SCVI.load(model_dir, adata=merged_hvg)

# ---- Obtain model outputs ----
SCVI_LATENT_KEY = "X_scVI"
latent = model.get_latent_representation()
merged_hvg.obsm[SCVI_LATENT_KEY] = latent

SCVI_NORMALIZED_KEY = "scvi_normalized"
merged_hvg.layers[SCVI_NORMALIZED_KEY] = model.get_normalized_expression(library_size=10e4)

# ---- Visualization without batch correction ----
sc.tl.pca(merged_hvg)
sc.pp.neighbors(merged_hvg, n_pcs=50, n_neighbors=20)
sc.tl.umap(merged_hvg, min_dist=0.3)
sc.pl.umap(merged_hvg, color=["sample"], frameon=False)

# ---- Visualization with batch correction (scVI latent space) ----
sc.pp.neighbors(merged_hvg, use_rep=SCVI_LATENT_KEY)
sc.tl.umap(merged_hvg, min_dist=0.3)
sc.pl.umap(merged_hvg, color=["sample"], frameon=False)

# ============================================================
# III. Clustering on the scVI latent space
# ============================================================

SCVI_CLUSTERS_KEY = "leiden_scVI"

for res in [0.1, 0.2, 0.5, 1.0]:
    sc.tl.leiden(merged_hvg, key_added=f"{SCVI_CLUSTERS_KEY}_{res:3.1f}", resolution=res)

sc.pl.umap(
    merged_hvg,
    color=['leiden_scVI_0.1', 'leiden_scVI_0.2', 'leiden_scVI_0.5', 'leiden_scVI_1.0'],
    legend_loc='on data',
    frameon=False,
)

# ---- Marker gene panel used for cluster identity ----
marker_genes = {
    "KC_suprabasal": ["Krt1", "Krt10", 'Flg', 'Cnfn'],
    "KC_basal": ["Krt5", "Krt14", "Col17a1"],
    "Fibroblasts": ["Col1a1", "Lum", "Pdgfra"],
    "Mast_cell": ['Fcer1a', 'Cd55', 'Gata2', 'Mrgprb1'],
    "Melanocyte": ["Mlana", "Tyrp1", "Dct"],
    "Nerve": ['Dag1', 'Cryab', 'Pmp22'],
    "Smooth_muscle": ["Tinagl1", "Mylk", "Tagln"],
    "Eccrine_gland": ["Mucl1", 'Ceacam5'],
    "Endo_vascular": ["Cd34", "Pecam1", "Vwf"],
    "Endo_lymphatic": ["Lyve1"],
    "T cells": ['Cd3d', "Trac", "Trbc1"],
    "CD4+ T": ["Cd4", "Il7r", "Cxcr4"],  # Cxcr4: T helper cells
    "CD8+ T": ["Cd8a"],
    "Th17": ['Rorc', 'Il17a', 'Il17f', 'Ccr6', 'Il23r', 'Stat3'],
    "NK": ['Nkg7', 'Gzmb', 'Gzmk', "Klrd1"],
    "DC": ['H2-Aa', 'H2-Ab1'],
    "Langerhans": ["Cd207", 'Epcam'],
    "Macrophage": ['Csf1r', 'Cd68', 'Mertk', 'Marco', 'Cd163', 'Adgre1'],
    "Neutrophil": ['Ly6g', 'Mpo', 'Elane', 'Ltf', 'Ngp'],
    "pDC": ['Bst2', 'Siglech', 'Tcf4'],
    'cDC1': ['Itgae', 'Xcr1', 'Cd8a'],
    "cDC2": ["Itgam", 'Sirpa', 'Cd4', "Cx3cr1"],
    "B cells": ['Cd79a'],
    "Plasma cells": ['Mzb1', 'Igkc', 'Ighg1'],
}

sc.pl.umap(
    merged_hvg,
    color=['Cd207', 'Epcam', 'H2-Aa', 'H2-Ab1', 'Cd74', 'Cd68',
           'Csf1r', 'Cd163', 'Cd3d', 'Cd3g', 'Il17a', 'Il22'],
)

sc.pl.dotplot(merged_hvg, marker_genes, groupby='leiden_scVI_0.2', standard_scale="var")


def extract_degs(adata, groupby, output_path, logfc_thresh=1.5, pval_thresh=0.05,
                  numeric_sort=True):
    """Run rank_genes_groups (Wilcoxon) and save filtered DEGs per cluster to CSV."""
    sc.tl.rank_genes_groups(adata, groupby=groupby, method="wilcoxon")
    rank_genes_groups = adata.uns['rank_genes_groups']
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

extract_degs(merged_hvg, 'leiden_scVI_0.2', 'results/cluster_markers_DEG/CD14_cluster_markers_DEG_0.2.csv')
extract_degs(merged_hvg, 'leiden_scVI_1.0', 'results/cluster_markers_DEG/CD14_cluster_markers_DEG_1.0.csv')

# ---- Sub-cluster the myeloid compartment (cluster 0 at resolution 0.2) at higher resolution ----
sc.tl.leiden(merged_hvg, resolution=0.4, restrict_to=['leiden_scVI_0.2', ['0']], key_added='myel')
sc.pl.umap(merged_hvg, color='myel', legend_loc='on data', frameon=False)

sc.pl.umap(merged, color=['Itgax', 'Adgre1'])
sc.pl.umap(merged, color=['Tcrg-V5', 'Adgre1'])

extract_degs(merged_hvg, 'myel', 'results/cluster_markers_DEG/CD14_cluster_markers_myel.csv',
             numeric_sort=False)

# ============================================================
# IV. Marker-based cluster annotation
# ============================================================

merged_hvg.obs["annot"] = merged_hvg.obs["myel"].map({
    '3': 'FB_1',
    '4': 'FB_2',
    '5': 'FB_3',

    '12': 'SM',
    '13': 'Endo_L',
    '7': 'Endo_V',
    '10': 'KC',
    '9': 'Melanocyte',

    '2': 'Neutrophil',

    '14': 'MastCell',
    '11': 'LC',
    '16': 'B',

    '1': 'T_IL17',
    '6': 'Th2',  # Gata3, Il13
    '15': 'CTL',

    '8': 'NKT',
    '0,2': 'DC',
    '0,4': 'DC',
    '0,5': 'DC',
    '0,6': 'DC',
    '0,7': 'DC',
    '0,8': 'DC',
    '0,9': 'DC',

    '0,0': 'Macrophage',
    '0,1': 'Macrophage',
    '0,3': 'Macrophage',
})

merged_hvg.obs['annot'] = np.where(
    merged_hvg.obs['leiden_scVI_0.5'] == '18', 'Nerve', merged_hvg.obs['annot']
)

merged_hvg.write('data/mouse_skin_merged_hvg.h5ad')

# ---- Copy lv0 annotation and UMAP embedding back onto the un-subsetted object ----
merged = merged[merged_hvg.obs_names].copy()

assert all(merged.obs_names == merged_hvg.obs_names), "Error: cell indices do not match!"

merged.obs['annot'] = merged_hvg.obs['annot'].copy()

for key in ['X_pca', 'X_umap', 'X_scVI']:
    merged.obsm[key] = merged_hvg.obsm[key].copy()

merged.write('data/mouse_skin_annotated.h5ad')

merged = sc.read('data/mouse_skin_annotated.h5ad')

sc.pl.umap(merged, color='annot', legend_loc='on data', frameon=False, legend_fontsize=9)
sc.pl.umap(merged, color='annot', legend_loc=None, frameon=False)

# ---- Non-overlapping cluster labels on UMAP ----
from adjustText import adjust_text

umap_coords = merged.obsm['X_umap']
merged.obs['umap1'] = umap_coords[:, 0]
merged.obs['umap2'] = umap_coords[:, 1]

sc.pl.umap(merged, color='annot', legend_loc=None, frameon=False, show=False)

cluster_centroids = merged.obs.groupby('annot')[['umap1', 'umap2']].mean()
texts = [plt.text(x, y, str(cluster), fontsize=9) for cluster, (x, y) in cluster_centroids.iterrows()]
adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
plt.show()

sc.pl.umap(merged, color='Cd14')

# ---- Final annotated dotplot (ordered by broad cell-type grouping) ----
final_marker_genes = {
    "DC": ['H2-Aa', 'H2-Ab1', 'Itgax'],
    "LC": ["Cd207", 'Epcam'],
    "Macrophage": ['Csf1r', 'Cd68', 'Adgre1'],
    "Neutrophil": ['Ly6g', 'Csf3r', 'Hcar2'],
    "T_IL17": ['Cd3e', 'Rorc', 'Il17a', 'Il17f', 'Ccr6', 'Il23r'],
    'CTL': ['Ccl5', 'Cd8a', 'Gzma'],
    'Th2': ['Gata3', 'Il13'],
    "NK": ['Nkg7', 'Gzmb', "Klrd1"],
    "B cells": ['Cd19', 'Iglc3', 'Mzb1', 'Cd79a'],

    "Mast_cell": ['Fcer1a', 'Cd55', 'Mrgprb1'],
    "Keratinocytes": ['Krt5', 'Krt14', "Krt1", "Krt10"],
    "Melanocyte": ["Mlana", "Tyrp1", "Dct"],
    "Nerve": ['Kcna1', 'Mpz', 'Cadm4'],
    "Smooth_muscle": ["Tinagl1", "Mylk", "Tagln"],
    "Endo_lymphatic": ["Lyve1"],
    "Endo_vascular": ["Cd34", "Pecam1", "Vwf"],

    "Fibroblasts": ['Lum', "Col3a1", 'Col23a1', 'Col6a1'],
}

desired_order = [
    'DC', 'LC', 'Macrophage', 'Neutrophil', 'T_IL17', 'CTL', 'Th2', 'NKT', 'B', 'MastCell',
    'KC', 'Melanocyte', 'Nerve', 'SM', 'Endo_L', 'Endo_V', 'FB_1', 'FB_2', 'FB_3',
]

merged.obs['annot'] = pd.Categorical(merged.obs['annot'], categories=desired_order, ordered=True)

sc.pl.dotplot(merged, final_marker_genes, groupby='annot', standard_scale="var")

# ============================================================
# V. Cluster cell count per group / sample
# ============================================================

df = merged.obs[['group', 'sample', 'annot']].copy()

cell_count_group = pd.crosstab(index=df['annot'], columns=df['group'])
print(cell_count_group)

cell_count_sample = pd.crosstab(index=df['annot'], columns=df['sample'])
print(cell_count_sample)
