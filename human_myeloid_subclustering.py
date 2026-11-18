"""
Human skin scRNA-seq: Myeloid subclustering

Starting from the level-0 annotated object (see human_integration_clustering.py),
subset to the 'Myeloid' compartment and re-integrate/re-cluster at higher resolution
using Harmony (instead of scVI), then annotate myeloid subtypes.

Pipeline overview:
  I.   Subset to Myeloid cells, HVG selection (with manual marker rescue)
  II.  Integration with Harmony + Leiden clustering at multiple resolutions
  III. Marker-based annotation (lv1) + DEG extraction
  IV.  Finalize annotated object (drop intermediate obs/uns columns)

Note: Data paths below are relative placeholders (data/, results/).
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

import matplotlib.pyplot as plt
import seaborn as sns

# ---- scanpy settings ----
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=75, color_map='OrRd')
sc.logging.print_version_and_date()

# ============================================================
# I. Subset to Myeloid cells, HVG selection
# ============================================================

merged = sc.read('data/human_skin_annotated.h5ad')
myel = merged[merged.obs['lv0'] == 'Myeloid'].copy()

sc.pp.highly_variable_genes(myel, min_mean=0.01, max_mean=5, min_disp=0.2)

hvg_list = sorted(myel.var[myel.var['highly_variable'] == True].index)
print(len(hvg_list))

# ---- Check specific myeloid markers are included in the HVGs ----
genes_to_check = [
    'CD207', 'CD1A', 'HLA-DQB2',
    'LAMP3', 'FSCN1',
    'ITGAX',
    'FCER1A', 'CLEC10A', 'CD1C', 'CD1E',
    'CLEC9A', 'WDFY4', 'CADM1', 'XCR1',
    'CLEC4C', 'IL3RA',
    'CD68',
    'CD14', 'CD209', 'ITGAM',
    'CD80',
    'FCN1', 'FCGR3A',
    'SPP1', 'MRC1', 'MS4A4A', 'CD163',
    'CSF3R', 'MPO', 'CD33',
]

genes_absent = [gene for gene in genes_to_check if gene not in hvg_list]
print(genes_absent)

# ---- Manually add omitted markers to HVGs ----
for gene in genes_absent:
    if gene in myel.var_names:
        myel.var.loc[gene, 'highly_variable'] = True

updated_hvg_genes = myel.var['highly_variable']
print(f"Total number of HVGs after appending: {updated_hvg_genes.sum()}")

# ============================================================
# II. Integration: Harmony
# ============================================================

sc.tl.pca(myel, use_highly_variable=True)
sc.pl.pca_variance_ratio(myel, n_pcs=50, log=True)

# Step 1: Harmony integration using scanpy's wrapper function
sce.pp.harmony_integrate(myel, key='SerSam', max_iter_harmony=50)

# Step 2: Neighbor graph construction using Harmony-corrected PCA embeddings
sc.pp.neighbors(myel, use_rep='X_pca_harmony', n_neighbors=15, n_pcs=50)

# Step 3: UMAP calculation
sc.tl.umap(myel)

palette = sns.color_palette('tab20', n_colors=myel.obs['SerSam'].nunique())
sc.pl.umap(myel, color='SerSam', frameon=False, palette=palette, legend_loc='none')

# ---- Marker gene panels used to identify myeloid subtypes ----
sc.pl.umap(myel, color=['CD207', 'CD1A', 'HLA-DQB2', 'ITGAX'])  # ITGAX: CD11c
sc.pl.umap(myel, color=['CLEC9A', 'WDFY4', 'CADM1', 'XCR1'])  # cDC1
sc.pl.umap(myel, color=['FCER1A', 'CLEC10A', 'CD1C', 'CD1E'])  # cDC2
sc.pl.umap(myel, color=['LAMP3', 'FSCN1', 'ICAM1', 'CCR7', 'CCL22', 'CCL19', 'IDO1'])  # matureDC/mregDC
sc.pl.umap(myel, color=['CLEC4C', 'IL3RA', 'IGKC', 'JCHAIN'])  # pDC
sc.pl.umap(myel, color=['CD68', 'FCN1', 'FCGR1A', 'CD40'])
sc.pl.umap(myel, color=['SPP1', 'MRC1', 'MS4A4A', 'CD163', 'IL10', 'TGFB1'])
sc.pl.umap(myel, color=['CSF3R', 'S100A8'], size=5)

# moDC
sc.pl.umap(myel, color=['CD14', 'CD209', 'FCGR3A', 'MAFB', 'CD83', 'CD1C'])

# DC3
sc.pl.umap(myel, color=['CD14', 'ITGAX', 'CD1C', 'CLEC10A', 'CD163', 'EREG', 'SLC2A3', 'IL23A'])

# ---- Leiden clustering at multiple resolutions ----
for res in [0.2, 0.5, 1.0, 2.0]:
    sc.tl.leiden(myel, key_added=f"leiden_res_{res:3.1f}", resolution=res)

sc.pl.umap(
    myel,
    color=['leiden_res_0.2', "leiden_res_0.5", "leiden_res_1.0", 'leiden_res_2.0'],
    frameon=False,
    legend_loc="on data",
)

for res in [3.0]:
    sc.tl.leiden(myel, key_added=f"leiden_res_{res:3.1f}", resolution=res)

sc.pl.umap(myel, color='leiden_res_3.0', legend_loc='on data', frameon=False)

myel.write('data/human_skin_myeloid.h5ad')

# ============================================================
# III. Annotation
# ============================================================


def extract_degs(adata, groupby, key_added, groups=None, reference='rest',
                  output_path=None, logfc_thresh=1.5, pval_thresh=0.05, numeric_sort=True):
    """Run rank_genes_groups (Wilcoxon) and optionally save filtered DEGs to CSV.

    If `groups` is given, only those group(s) are compared against `reference`
    (e.g. one cluster vs. another, to help resolve an ambiguous split).
    Otherwise all groups in `groupby` are compared against the rest.
    """
    sc.tl.rank_genes_groups(
        adata, groupby=groupby, groups=groups if groups else 'all',
        reference=reference, key_added=key_added, method="wilcoxon",
    )
    rank_genes_groups = adata.uns[key_added]
    result_groups = groups if groups else rank_genes_groups['names'].dtype.names
    if numeric_sort and not groups:
        result_groups = sorted(result_groups, key=lambda x: int(x))

    result_dfs = []
    for group in result_groups:
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
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        all_degs.to_csv(output_path, index=False)
    return all_degs


os.makedirs('results/cluster_markers_myeloid', exist_ok=True)

# ---- Cluster-specific DE genes at resolution 0.2 (used for annotation below) ----
extract_degs(
    myel, groupby="leiden_res_0.2", key_added="rank_genes_leiden_0.2",
    output_path="results/cluster_markers_myeloid/cluster_markers_lv1_myel_deg_0.2.csv",
)

# ---- Pairwise comparison to resolve an ambiguous split (cluster 22 vs. 15 at res 2.0) ----
extract_degs(
    myel, groupby="leiden_res_2.0", groups=['22'], reference='15', key_added="sub",
    output_path="results/cluster_markers_myeloid/DEG_cluster22_vs_cluster15_res2.0.csv",
)

sc.pl.umap(myel, color='pct_counts_mt')

# ---- lv1 annotation ----
myel.obs["lv1_v401_nomt"] = myel.obs["leiden_res_0.2"].map({
    '0': 'cDC2',       # CLEC10A, FCER1A, CST7, CD1C
    '1': 'M1mac',      # CXCL2, CXCL3, CXCL8, EREG
    '2': 'M2mac',      # C1QA-C, MS4A6A, CD68, CD14, MAFB, CD163, MS4A4A
    '3': 'cDC1',       # CLEC9A, WDFY4, CADM1, XCR1
    '4': 'mregDC',     # LAMP3, FSCN1, CCR7, CCL22, IDO1
    '5': 'pDC',        # IGKC, JCHAIN, GZMB, IL3RA, MZB1
    '6': 'LC',         # HLA-DQB2, CD207, CD1A
})

myel.obs['lv1_v401_nomt'] = np.where(
    myel.obs['leiden_res_0.5'] == '1', 'M2mac', myel.obs['lv1_v401_nomt']
)
myel.obs['lv1_v401_nomt'] = np.where(
    myel.obs['leiden_res_1.0'] == '12', 'Neutrophil', myel.obs['lv1_v401_nomt']
)
myel.obs['lv1_v401_nomt'] = np.where(
    myel.obs['leiden_res_2.0'] == '15', 'moDC', myel.obs['lv1_v401_nomt']
)
myel.obs['lv1_v401_nomt'] = np.where(
    myel.obs['leiden_res_2.0'] == '17', 'mregDC', myel.obs['lv1_v401_nomt']
)

sc.pl.umap(myel, color='lv1_v401_nomt', ncols=1, size=5)

myel.write('data/human_skin_myeloid.h5ad')

marker_genes = {
    'LC': ['CD207', 'CD1A', 'HLA-DQB2'],
    'M1mac': ['FCN1'],
    'M2mac': ['CD163', 'MRC1', 'MS4A6A', 'C1QA'],
    'Neutrophil': ['CSF3R', 'S100A8', 'S100A9'],
    'cDC1': ['CLEC9A', 'XCR1', 'WDFY4', 'CADM1'],
    'cDC2': ['CLEC10A', 'FCER1A', 'CD1C'],
    'moDC': ['CD14', 'FCGR3A'],
    'mregDC': ['FSCN1', 'LAMP3', 'CCL22', 'IDO1'],
    'pDC': ['IGKC', 'JCHAIN', 'IL3RA', 'CLEC4C'],
}

sc.pl.dotplot(myel, marker_genes, groupby='lv1_v401_nomt', standard_scale="var")

# ============================================================
# IV. Finalize annotated object
# ============================================================

myel_final = myel.copy()
myel_final.obs['lv1'] = myel_final.obs['lv1_v401_nomt'].copy()

# Drop intermediate per-resolution obs columns, keeping only the final 'lv1' label
obs_cols_to_drop = [
    'Unnamed: 0', 'leiden_res_0.2', 'leiden_res_0.5',
    'leiden_res_1.0', 'leiden_res_2.0', 'leiden_res_2.5', 'leiden_res_3.0',
]
for col in obs_cols_to_drop:
    if col in myel_final.obs:
        del myel_final.obs[col]

# Drop intermediate uns entries (clustering colors, per-resolution DE results, etc.)
uns_keys_to_drop = [
    'c11', 'cdm_leiden_res_0.2', 'cdm_leiden_res_0.5',
    'leiden_res_0.2', 'leiden_res_0.2_colors', 'leiden_res_0.5', 'leiden_res_0.5_colors',
    'leiden_res_1.0', 'leiden_res_1.0_colors', 'leiden_res_2.0', 'leiden_res_2.0_colors',
    'leiden_res_2.5', 'leiden_res_2.5_colors', 'leiden_res_3.0', 'leiden_res_3.0_colors',
    'predicted_doublet_colors', 'rank_genes_leiden_0.2', 'rank_genes_leiden_0.5',
    'rank_genes_leiden_1.0', 'rank_genes_leiden_2.0', 'rank_genes_leiden_3.0', 'sub',
]
for key in uns_keys_to_drop:
    if key in myel_final.uns:
        del myel_final.uns[key]

myel_final.write('data/human_skin_myeloid_final.h5ad')
