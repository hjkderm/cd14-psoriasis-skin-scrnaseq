"""
Batch-effect evaluation for the integrated human skin datasets

Assesses how well integration removed study-level batch effects, at two levels:
  Part 1. Whole tissue  — scVI integration    (see human_integration_clustering.py)
  Part 2. Myeloid subset — Harmony integration (see human_myeloid_subclustering.py)

For each level:
  - UMAP before integration (merge-only PCA) vs. after integration, colored by
    the batch variable (Series = independent GEO dataset; SerSam = sample)
  - iLISI (batch mixing) and cLISI (cell-type separation) via the scib package,
    using dataset (Series) as the batch variable

Note: Data paths below are relative placeholders (data/).
Adjust them to match your own directory layout before running.
"""

from warnings import filterwarnings
filterwarnings('ignore')

import numpy as np
import pandas as pd

import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

import scib

BATCH_KEY = 'Series'   # 23 independent GEO datasets

# ============================================================
# Part 1. Whole tissue — scVI
# ============================================================

pso = sc.read('data/human_skin_annotated.h5ad')

# ---- Compute the "before integration" embedding (merge-only PCA) ----
# The integrated (scVI) UMAP already lives in .obsm['X_umap']; back it up before
# overwriting, then recompute a UMAP from the uncorrected PCA for comparison.
if 'X_umap' in pso.obsm:
    pso.obsm['X_umap_scvi'] = pso.obsm['X_umap'].copy()

sc.pp.neighbors(pso, n_pcs=50, n_neighbors=15, key_added='neighbors_merge_only')
sc.tl.umap(pso, min_dist=0.3, neighbors_key='neighbors_merge_only')
pso.obsm['X_umap_merge_only'] = pso.obsm['X_umap'].copy()

# Restore the scVI UMAP as the default embedding
pso.obsm['X_umap'] = pso.obsm['X_umap_scvi'].copy()

# ---- Batch-colored UMAPs, before vs. after ----
palette = sns.color_palette('tab20', n_colors=pso.obs['SerSam'].nunique())

sc.pl.embedding(pso, basis='X_umap_merge_only', color='SerSam',
                title='Merge only (no scVI)', palette=palette,
                frameon=False, legend_loc='none')
sc.pl.embedding(pso, basis='X_umap_merge_only', color=BATCH_KEY,
                title='Merge only (no scVI)', frameon=False)

sc.pl.embedding(pso, basis='X_umap_scvi', color='SerSam',
                title='After scVI', palette=palette,
                frameon=False, legend_loc='none')
sc.pl.embedding(pso, basis='X_umap_scvi', color=BATCH_KEY,
                title='After scVI', frameon=False)

pso.write('data/human_skin_annotated.h5ad')

# ---- Restore the scVI latent representation if it is not on the object ----
# iLISI/cLISI are computed on the latent space, not the UMAP. If X_scVI is
# missing (e.g. it was dropped when the object was subset), copy it back from
# the object the scVI model was trained on, aligned by cell barcode.
if 'X_scVI' not in pso.obsm:
    merged = sc.read('data/human_skin_merged_hvg.h5ad')

    common_cells = pso.obs_names.intersection(merged.obs_names)
    print(f"Matching cells: {len(common_cells)} / {pso.n_obs}")

    latent_df = pd.DataFrame(merged.obsm['X_scVI'], index=merged.obs_names)
    latent_aligned = latent_df.reindex(pso.obs_names)

    n_missing = latent_aligned.isna().any(axis=1).sum()
    print(f"Cells without matching scVI latent: {n_missing}")

    pso.obsm['X_scVI'] = latent_aligned.values
    print("X_scVI restored:", pso.obsm['X_scVI'].shape)

    del merged
    import gc
    gc.collect()

    pso.write('data/human_skin_annotated.h5ad')

# ---- iLISI / cLISI (whole tissue) ----
# iLISI: higher = better mixing of batches. cLISI: closer to 1 = better
# preservation of cell-type separation.
label_key = 'lv0'

ilisi_before = scib.metrics.ilisi_graph(
    pso, batch_key=BATCH_KEY, type_='embed', use_rep='X_pca', k0=90
)
ilisi_after = scib.metrics.ilisi_graph(
    pso, batch_key=BATCH_KEY, type_='embed', use_rep='X_scVI', k0=90
)
print(f"iLISI before (merge-only PCA): {ilisi_before:.3f}")
print(f"iLISI after  (scVI):           {ilisi_after:.3f}")

clisi_before = scib.metrics.clisi_graph(
    pso, label_key=label_key, type_='embed', use_rep='X_pca', k0=90
)
clisi_after = scib.metrics.clisi_graph(
    pso, label_key=label_key, type_='embed', use_rep='X_scVI', k0=90
)
print(f"cLISI before (merge-only PCA): {clisi_before:.3f}")
print(f"cLISI after  (scVI):           {clisi_after:.3f}")

# ============================================================
# Part 2. Myeloid subset — Harmony
# ============================================================

myel = sc.read('data/human_skin_myeloid.h5ad')

# ---- Keep Series colors consistent with the whole-tissue plots ----
categories_pso = pso.obs[BATCH_KEY].cat.categories.tolist()
color_map = dict(zip(categories_pso, pso.uns[f'{BATCH_KEY}_colors']))

myel.obs[BATCH_KEY] = myel.obs[BATCH_KEY].astype('category')
present = [c for c in categories_pso if c in myel.obs[BATCH_KEY].unique()]
myel.obs[BATCH_KEY] = myel.obs[BATCH_KEY].cat.set_categories(present)
myel.uns[f'{BATCH_KEY}_colors'] = [color_map[c] for c in present]

# ---- Compute the "before integration" embedding (merge-only PCA) ----
sc.pp.neighbors(myel, use_rep='X_pca', key_added='neighbors_merge_only')

if 'X_umap' in myel.obsm:
    myel.obsm['X_umap_harmony'] = myel.obsm['X_umap'].copy()

sc.tl.umap(myel, neighbors_key='neighbors_merge_only')
myel.obsm['X_umap_merge_only'] = myel.obsm['X_umap'].copy()

# Restore the Harmony UMAP as the default embedding
myel.obsm['X_umap'] = myel.obsm['X_umap_harmony'].copy()

myel.write('data/human_skin_myeloid.h5ad')

# ---- Check for overplotting: cell counts are very uneven across datasets ----
series_counts = myel.obs.groupby(BATCH_KEY).size().sort_values(ascending=False)
print(series_counts)

top_series = series_counts.index[0]
print(f"Top series: {top_series}")
print(myel.obs.loc[myel.obs[BATCH_KEY] == top_series,
                   ['Lesion', 'Disease', 'BodySite', 'Category']].drop_duplicates())

# Is the largest dataset skewed towards a particular myeloid subtype?
target_celltypes = ['cDC1', 'cDC2', 'LC', 'moDC', 'M1mac', 'M2mac']
myel_obs_filtered = myel.obs[myel.obs['lv1'].isin(target_celltypes)].copy()
myel_obs_filtered['lv1'] = myel_obs_filtered['lv1'].cat.remove_unused_categories()

dataset_prop = myel_obs_filtered.groupby([BATCH_KEY, 'lv1']).size().unstack(fill_value=0)
dataset_prop = dataset_prop.div(dataset_prop.sum(axis=1), axis=0)
print(dataset_prop)

# ---- Batch-colored UMAPs, before vs. after ----
# Cells per dataset are capped for VISUALIZATION ONLY, so that a few very large
# datasets do not visually dominate the plot. All metrics above and all
# statistical analyses elsewhere use the full data.
np.random.seed(0)
max_per_series = 1000

sampled_indices = []
for series in myel.obs[BATCH_KEY].unique():
    idx = np.where(myel.obs[BATCH_KEY] == series)[0]
    if len(idx) > max_per_series:
        idx = np.random.choice(idx, max_per_series, replace=False)
    sampled_indices.extend(idx)

myel_balanced_plot = myel[sampled_indices].copy()

# Shuffle so that no dataset is systematically drawn on top of the others
shuffle_idx = np.random.permutation(myel_balanced_plot.n_obs)
myel_balanced_plot = myel_balanced_plot[shuffle_idx].copy()

sc.pl.embedding(myel_balanced_plot, basis='X_umap_merge_only', color=BATCH_KEY,
                title=f'Merge only (no Harmony), capped at {max_per_series} cells/dataset',
                size=8, alpha=0.6, legend_loc='right margin', frameon=False)
sc.pl.embedding(myel_balanced_plot, basis='X_umap_harmony', color=BATCH_KEY,
                title=f'After Harmony, capped at {max_per_series} cells/dataset',
                size=8, alpha=0.6, legend_loc='right margin', frameon=False)

# ---- iLISI / cLISI (myeloid) ----
label_key = 'lv1'

ilisi_before = scib.metrics.ilisi_graph(
    myel, batch_key=BATCH_KEY, type_='embed', use_rep='X_pca', k0=90
)
ilisi_after = scib.metrics.ilisi_graph(
    myel, batch_key=BATCH_KEY, type_='embed', use_rep='X_pca_harmony', k0=90
)
print(f"iLISI before (merge-only PCA): {ilisi_before:.3f}")
print(f"iLISI after  (Harmony):        {ilisi_after:.3f}")

clisi_before = scib.metrics.clisi_graph(
    myel, label_key=label_key, type_='embed', use_rep='X_pca', k0=90
)
clisi_after = scib.metrics.clisi_graph(
    myel, label_key=label_key, type_='embed', use_rep='X_pca_harmony', k0=90
)
print(f"cLISI before (merge-only PCA): {clisi_before:.3f}")
print(f"cLISI after  (Harmony):        {clisi_after:.3f}")
