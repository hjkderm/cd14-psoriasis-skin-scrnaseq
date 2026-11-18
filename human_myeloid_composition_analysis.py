"""
Human myeloid compartment: composition robustness and marker expression checks

Part 1. Sample-size downsampling
        The public cohort is imbalanced (83 healthy vs 18 psoriasis samples).
        To test whether myeloid subtype proportions are driven by that
        imbalance, healthy samples are repeatedly downsampled to match the
        psoriasis sample number and the proportions are recalculated.

Part 2. Per-dataset composition breakdown
        Myeloid subtype composition per GEO dataset, to check whether the
        overall composition is driven by any single dataset.

Part 3. FAS and CFLAR expression
        Mean expression and fraction of expressing cells per myeloid subtype,
        for all samples and separately for healthy and psoriasis samples. The
        fraction is reported alongside the mean, because mean expression alone
        is vulnerable to dropout for sparsely detected genes.

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

pso = sc.read('data/human_skin_annotated.h5ad')

target_celltypes = ['cDC1', 'cDC2', 'LC', 'moDC', 'M1mac', 'M2mac']

# Boolean mask over the full object; avoids subsetting/copying the whole AnnData
myeloid_mask = pso.obs['lv0'] == 'Myeloid'

# ============================================================
# Part 1. Sample-level downsampling of the healthy cohort
# ============================================================

# ---- Sample-level metadata ----
sample_meta = pso.obs[['SerSam', 'Series', 'Disease']].drop_duplicates().set_index('SerSam')
healthy_samples_all = sample_meta.loc[sample_meta['Disease'] == 'HC'].index.unique()
pso_samples_all = sample_meta.loc[sample_meta['Disease'] == 'Pso'].index.unique()
n_target = len(pso_samples_all)
print(f"Healthy samples: {len(healthy_samples_all)}, Psoriasis samples: {len(pso_samples_all)}")

# ---- Keep only healthy samples with enough myeloid cells to contribute ----
myeloid_counts = pso.obs.loc[myeloid_mask].groupby('SerSam').size()
min_myeloid_cells = 10
valid_healthy = [s for s in healthy_samples_all if myeloid_counts.get(s, 0) >= min_myeloid_cells]
print(f"Valid healthy samples after myeloid count filter: {len(valid_healthy)} / {len(healthy_samples_all)}")

# ---- Repeatedly downsample healthy to match psoriasis, recompute proportions ----
myeloid_obs = pso.obs.loc[myeloid_mask, ['SerSam', 'Series', 'lv1', 'Disease']]
myeloid_obs = myeloid_obs[myeloid_obs['lv1'].isin(target_celltypes)]

np.random.seed(0)
n_iterations = 50
results = []

for i in range(n_iterations):
    selected_healthy = np.random.choice(valid_healthy, size=n_target, replace=False)
    keep_samples = list(selected_healthy) + list(pso_samples_all)

    sub_obs = myeloid_obs[myeloid_obs['SerSam'].isin(keep_samples)]

    prop = sub_obs.groupby(['Disease', 'lv1']).size().unstack(fill_value=0)
    prop = prop.div(prop.sum(axis=1), axis=0)
    prop['iteration'] = i
    results.append(prop)

all_results = pd.concat(results).reset_index()

# ---- Reference: proportions from the full (83 vs 18) cohort ----
orig_prop = myeloid_obs.groupby(['Disease', 'lv1']).size().unstack(fill_value=0)
orig_prop = orig_prop.div(orig_prop.sum(axis=1), axis=0)
print("=== Original (83 HC vs 18 Pso) proportion ===")
print(orig_prop)

# ---- Downsampled distribution vs. the full-cohort reference ----
melted = all_results.melt(id_vars=['Disease', 'iteration'], value_vars=target_celltypes,
                          var_name='lv1', value_name='proportion')

fig, ax = plt.subplots(figsize=(8, 5))
sns.boxplot(data=melted, x='lv1', y='proportion', hue='Disease', ax=ax)

# Overlay the full-cohort proportions as diamonds, as a comparison baseline
for j, ct in enumerate(target_celltypes):
    for disease, color in zip(['HC', 'Pso'], ['blue', 'orange']):
        if disease in orig_prop.index and ct in orig_prop.columns:
            ax.scatter(j, orig_prop.loc[disease, ct], marker='D', color=color,
                       edgecolor='black', s=60, zorder=5)

ax.set_title(f'Myeloid subtype proportion across {n_iterations} downsampling iterations\n'
             f'(Healthy n={n_target} matched to Psoriasis; diamonds = original 83 vs 18 result)')
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

# ============================================================
# Part 2. Per-dataset composition breakdown
# ============================================================
# Checks whether the overall composition is a batch-driven artifact of any
# single dataset.

myeloid_obs['lv1'] = myeloid_obs['lv1'].cat.remove_unused_categories()

dataset_prop = myeloid_obs.groupby(['Series', 'lv1']).size().unstack(fill_value=0)
dataset_prop = dataset_prop.div(dataset_prop.sum(axis=1), axis=0)

fig, ax = plt.subplots(figsize=(10, 6))
dataset_prop.plot(kind='bar', stacked=True, ax=ax, colormap='tab10')
ax.set_ylabel('Proportion')
ax.set_title('Myeloid subtype composition per dataset (Series)')
ax.legend(title='lv1', bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.xticks(rotation=90)
fig.subplots_adjust(right=0.78)
plt.show()

# ---- Number of samples per dataset, split by disease status ----
series_sample_counts = (
    pso.obs.loc[myeloid_mask, ['SerSam', 'Series', 'Disease']]
    .drop_duplicates()
    .groupby(['Series', 'Disease']).size().unstack(fill_value=0)
)
print("=== Number of samples per Series, by Disease ===")
print(series_sample_counts)

# ---- Same breakdown, split by disease, marking datasets with no myeloid cells ----
# Some datasets contribute samples but no myeloid cells at all; those are shown
# as hatched grey bars rather than silently omitted.
all_series_disease = pso.obs[['Series', 'Disease']].drop_duplicates()
all_series_disease = all_series_disease.set_index(['Series', 'Disease']).sort_index().index

dataset_prop_sd = myeloid_obs.groupby(['Series', 'Disease', 'lv1']).size().unstack(fill_value=0)
dataset_prop_sd = dataset_prop_sd.div(dataset_prop_sd.sum(axis=1), axis=0)

dataset_prop_full = dataset_prop_sd.reindex(all_series_disease)

no_myeloid_combos = dataset_prop_full[dataset_prop_full.isna().all(axis=1)].index
print("Series/Disease combinations with samples but no myeloid cells:")
print(list(no_myeloid_combos))

dataset_prop_plot = dataset_prop_full.fillna(0)

fig, ax = plt.subplots(figsize=(14, 6))
dataset_prop_plot.plot(kind='bar', stacked=True, ax=ax, colormap='tab10', width=0.8)

x_positions = {idx: i for i, idx in enumerate(dataset_prop_plot.index)}
for combo in no_myeloid_combos:
    ax.bar(x_positions[combo], 1.0, color='lightgrey', hatch='//',
           edgecolor='black', width=0.8, zorder=0)

ax.set_xticklabels([f"{series}\n({disease})" for series, disease in dataset_prop_plot.index],
                   rotation=90, fontsize=8)
ax.set_ylabel('Proportion')
ax.set_title('Myeloid subtype composition per dataset (Series), split by Disease\n'
             '(hatched grey = sample present but no myeloid cells captured)')

handles, labels = ax.get_legend_handles_labels()
handles.append(plt.Rectangle((0, 0), 1, 1, facecolor='lightgrey', hatch='//', edgecolor='black'))
labels.append('No myeloid captured')
ax.legend(handles, labels, title='lv1', bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

fig.subplots_adjust(right=0.78, bottom=0.25)
plt.show()

# ============================================================
# Part 3. CFLAR expression proportion
# ============================================================
# Mean expression alone is vulnerable to dropout for sparsely detected genes,
# so the fraction of cells with non-zero counts is reported alongside it.

gene = 'CFLAR'

expr = pso[myeloid_mask, gene].X.toarray().flatten()

df = pd.DataFrame({
    'lv1': pso.obs.loc[myeloid_mask, 'lv1'].values,
    'Disease': pso.obs.loc[myeloid_mask, 'Disease'].values,
    'expressing': expr > 0,
    'expr': expr,
})
df = df[df['lv1'].isin(target_celltypes)]

summary = pd.DataFrame({
    'pct_expressing_all': df.groupby('lv1', observed=True)['expressing'].mean().loc[target_celltypes],
    'pct_expressing_HC': df[df['Disease'] == 'HC'].groupby('lv1', observed=True)['expressing'].mean().loc[target_celltypes],
    'pct_expressing_Pso': df[df['Disease'] == 'Pso'].groupby('lv1', observed=True)['expressing'].mean().loc[target_celltypes],
    'mean_expr_all': df.groupby('lv1', observed=True)['expr'].mean().loc[target_celltypes],
})
print(summary)

# ---- FAS and CFLAR by disease status (mean expression and % expressing cells) ----
for g in ['FAS', 'CFLAR']:
    e = pso[myeloid_mask, g].X.toarray().flatten()
    d = pd.DataFrame({'lv1': pso.obs.loc[myeloid_mask, 'lv1'].values,
                      'Disease': pso.obs.loc[myeloid_mask, 'Disease'].values,
                      'expr': e})
    d = d[d['lv1'].isin(target_celltypes)]
    out = d.groupby(['Disease', 'lv1'], observed=True)['expr'].agg(
        mean='mean', pct_expressing=lambda x: (x > 0).mean() * 100)
    print(f'=== {g} by disease ===')
    print(out.round(3))

# ---- Dotplot of the core molecules across myeloid subtypes ----
myeloid_sub = pso[myeloid_mask & pso.obs['lv1'].isin(target_celltypes)].copy()
myeloid_sub.obs['lv1'] = myeloid_sub.obs['lv1'].cat.remove_unused_categories()

sc.pl.dotplot(
    myeloid_sub,
    var_names=['CD14', 'FAS', 'FADD', 'CFLAR'],
    groupby='lv1',
    categories_order=target_celltypes,
    standard_scale='var',   # scale each gene to 0-1 for comparability
    save='_core_molecules_dotplot.pdf',
)
