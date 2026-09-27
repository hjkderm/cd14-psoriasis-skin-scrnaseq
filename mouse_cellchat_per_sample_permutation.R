## ---------------------------------------------------------------------------
## Per-sample CellChat analysis + exact permutation test
##
## Rationale:
##   The main CellChat script (mouse_cellchat.R) subsets the Seurat object by
##   GROUP (KO / WT), pooling 3 mice into a single CellChat object per genotype.
##   This collapses n=3 biological replicates into n=1 pooled dataset, so
##   CellChat's internal permutation test (which shuffles CELL labels within one
##   dataset) cannot tell whether a genotype difference is reproducible across
##   mice.
##
##   This script instead:
##     (1) Subsets by 'sample' (6 individual mice) -> 6 independent CellChat
##         objects, preserving true n=3 per genotype.
##     (2) Extracts one scalar "pathway strength" value per mouse per pathway.
##     (3) Runs an EXACT permutation test (all C(6,3)=20 group-label
##         reassignments) comparing KO vs WT per pathway, with BH correction
##         for multiple pathways tested.
##     (4) Reports a direction-concordance metric as supporting evidence.
##     (5) Builds the supplementary table and the per-pathway dot + bar figure.
## ---------------------------------------------------------------------------

library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)
library(future)
options(stringsAsFactors = FALSE)

dir.create("results", showWarnings = FALSE)

## ---------------------------------------------------------------------------
## Part 0. Load data (same converted Seurat object used by mouse_cellchat.R)
## ---------------------------------------------------------------------------
obj <- readRDS("data/mouse_skin_for_cellchat.rds")
obj$samples <- obj$sample  # keep naming consistent with CellChat's expected 'samples' field

# Fix cell-type ('annot') factor levels ONCE, across the whole object, so that
# every per-sample CellChat object uses the same set of cell-group labels.
# (Some rare cell types may have 0 cells in a given mouse -- that's fine, it
# just means that group contributes 0 to the pathway sum for that mouse.)
all_levels <- levels(factor(obj$annot))
obj$annot  <- factor(obj$annot, levels = all_levels)

CellChatDB     <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB)  # same DB choice as the per-group script

sample_ids <- sort(unique(obj$sample))
cat("Samples found:\n"); print(sample_ids)
# Sanity check: confirm 3 KO + 3 WT
print(table(unique(data.frame(sample = obj$sample, group = obj$group))$group))

options(future.globals.maxSize = 2000 * 1024^2)
future::plan("multisession", workers = 4)

## ---------------------------------------------------------------------------
## Part 1. Run CellChat independently for each of the 6 mice
## ---------------------------------------------------------------------------
cc_list <- list()

for (s in sample_ids) {

  message("=== Running CellChat for sample: ", s, " ===")

  sub <- subset(obj, subset = sample == s)

  # Same construction as the per-group script, just per-mouse now
  cc <- createCellChat(object = sub, group.by = "annot", assay = "originalexp")
  cc@DB <- CellChatDB.use

  cc <- subsetData(cc)                      # subset signaling genes
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)

  cc <- computeCommunProb(cc, type = "triMean")
  cc <- filterCommunication(cc, min.cells = 10)   # NOTE: with fewer cells per mouse,
                                                  # some rare cell-type pairs may drop out;
                                                  # this is expected and is reported.

  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)

  cc_list[[s]] <- cc
  saveRDS(cc, file = file.path("results", paste0("cellchat_perSample_", s, ".rds")))
}

saveRDS(cc_list, file = "results/cellchat_perSample_list.rds")

## ---------------------------------------------------------------------------
## Part 2. Extract one "pathway strength" scalar per mouse per pathway
## ---------------------------------------------------------------------------
# Definition: for each pathway, sum the communication probability (netP$prob)
# over all sender x receiver cell-group pairs. This mirrors CellChat's own
# notion of pathway-level "information flow" (as used in rankNet), just kept
# per-sample instead of pooled.

pathway_strength <- function(cc) {
  netP <- cc@netP$prob            # array [sender, receiver, pathway]
  if (is.null(netP)) return(NULL)
  apply(netP, 3, sum)             # named numeric vector: pathway -> summed prob
}

strength_list <- lapply(cc_list, pathway_strength)

all_pathways <- unique(unlist(lapply(strength_list, names)))

strength_mat <- sapply(strength_list, function(x) {
  v <- setNames(rep(0, length(all_pathways)), all_pathways)
  v[names(x)] <- x
  v
})
# strength_mat: rows = pathway, columns = sample (6 mice). A pathway with 0 in
# a given mouse means it was not detected/passed filtering in that mouse --
# keep it as 0, don't drop the row, since that's informative for reproducibility.

meta_sample <- data.frame(sample = colnames(strength_mat))
group_lookup <- unique(data.frame(sample = obj$sample, group = obj$group))
meta_sample$group <- group_lookup$group[match(meta_sample$sample, group_lookup$sample)]

write.csv(strength_mat, "results/pathway_strength_perSample.csv")
write.csv(meta_sample,  "results/sample_group_map.csv", row.names = FALSE)

rm(obj)
gc()

## ---------------------------------------------------------------------------
## Part 3. Exact permutation test, KO (n=3) vs WT (n=3), per pathway
## ---------------------------------------------------------------------------
# With only 6 mice total (3 vs 3), there are exactly choose(6,3) = 20 distinct
# ways to relabel which 3 mice are "KO". This means an EXACT (exhaustive)
# permutation test is feasible -- and is more defensible than a Monte Carlo
# approximation, since every possible relabeling can be enumerated.
#
# IMPORTANT CAVEAT: with only 20 possible relabelings, the smallest achievable
# two-sided p-value is 2/20 = 0.1 (or 1/20 = 0.05 for a one-sided test). This
# is a hard mathematical floor, not a bug -- with n=3 per group, no test
# (parametric or non-parametric) can produce p < 0.05 two-sided from group-label
# permutation alone. Conclusions are therefore framed around effect consistency
# (direction concordance across mice) rather than the p-value alone.

perm_test_exact <- function(x, group) {
  # x: named numeric vector (per-sample pathway strength)
  # group: named character/factor vector, same names as x, values "KO"/"WT"
  stopifnot(all(names(x) == names(group)))
  n     <- length(x)
  n_ko  <- sum(group == "KO")

  obs_diff <- mean(x[group == "KO"]) - mean(x[group == "WT"])

  idx_combos <- combn(n, n_ko)  # all ways to choose which n_ko samples are "KO"
  null_diffs <- apply(idx_combos, 2, function(ko_idx) {
    mean(x[ko_idx]) - mean(x[-ko_idx])
  })

  p_two_sided <- mean(abs(null_diffs) >= abs(obs_diff) - 1e-12)  # tiny tol for float compare

  list(obs_diff = obs_diff, p_value = p_two_sided, n_perms = ncol(idx_combos))
}

group_vec <- meta_sample$group
names(group_vec) <- meta_sample$sample

pathway_results <- lapply(rownames(strength_mat), function(p) {
  x   <- strength_mat[p, ]
  res <- perm_test_exact(x, group_vec[names(x)])
  data.frame(
    pathway            = p,
    mean_KO            = mean(x[group_vec[names(x)] == "KO"]),
    mean_WT            = mean(x[group_vec[names(x)] == "WT"]),
    obs_diff_KO_vs_WT  = res$obs_diff,
    p_value_exact_perm = res$p_value,
    n_perms            = res$n_perms
  )
})

pathway_results_df <- do.call(rbind, pathway_results)

# Multiple-comparison correction across ALL pathways tested (BH/FDR)
pathway_results_df$p_adj_BH <- p.adjust(pathway_results_df$p_value_exact_perm, method = "BH")

pathway_results_df <- pathway_results_df[order(pathway_results_df$p_value_exact_perm), ]

## ---------------------------------------------------------------------------
## Part 4. Reproducibility check per pathway (supporting evidence)
## ---------------------------------------------------------------------------
# Because n=3 vs 3 permutation p-values are floor-limited (see caveat above),
# also report a simple concordance metric: does the KO vs WT direction of
# effect agree across ALL pairwise KO-WT mouse comparisons (3x3 = 9 pairs)?
# This is useful supporting evidence even for pathways that don't reach
# p < 0.1 in the exact permutation test.

concordance_check <- function(x, group) {
  ko_vals <- x[group == "KO"]
  wt_vals <- x[group == "WT"]
  pairwise_diffs <- outer(ko_vals, wt_vals, "-")  # 3x3 matrix of all KO-WT diffs
  mean(sign(pairwise_diffs) == sign(mean(pairwise_diffs)))  # fraction agreeing in direction
}

pathway_results_df$direction_concordance <- sapply(pathway_results_df$pathway, function(p) {
  x <- strength_mat[p, ]
  concordance_check(x, group_vec[names(x)])
})

write.csv(pathway_results_df, "results/pathway_permutationTest_results.csv", row.names = FALSE)

cat("\nDone. Key outputs:\n",
    " - pathway_strength_perSample.csv      (per-mouse pathway strengths)\n",
    " - sample_group_map.csv                (sample -> group mapping)\n",
    " - pathway_permutationTest_results.csv (exact permutation test + BH-adjusted p + concordance)\n")

## ---------------------------------------------------------------------------
## Part 5. Supplementary table: pathways discussed in the text
## ---------------------------------------------------------------------------
# Pathways included: those named in the Results text (WT-enriched ECM pathways,
# Cd14-/--enriched migration/adhesion/chemokine pathways) plus the
# antigen-presentation/complement pathways that were REMOVED from the text --
# included here for full transparency about why they were removed.
pathways_for_table <- c(
  # WT-enriched (ECM remodeling)
  "LAMININ", "FN1", "COLLAGEN", "TENASCIN",
  # Cd14-/--enriched (migration / adhesion / chemokine)
  "CXCL", "CCL", "ICAM",
  # Removed from text -- shown here for transparency
  "MHC-I", "CD80", "CD86", "COMPLEMENT"
)

# Check all requested pathways actually exist in the results (catch typos early)
missing_pw <- setdiff(pathways_for_table, pathway_results_df$pathway)
if (length(missing_pw) > 0) {
  warning("These pathways were not found in pathway_results_df: ",
          paste(missing_pw, collapse = ", "))
}

supp_table <- pathway_results_df[pathway_results_df$pathway %in% pathways_for_table, ]

# Order rows to match the order they're mentioned in text, rather than by p-value
supp_table <- supp_table[match(pathways_for_table, supp_table$pathway), ]

# Round numeric columns for readability. Mean/Difference columns use 3
# significant figures rather than a fixed 2 decimal places, since several
# pathways (e.g., MHC-I, CD80) have small absolute values that would
# otherwise all round to 0.00 and look identical.
sig3 <- function(x) signif(x, 3)

supp_table_out <- data.frame(
  Pathway                     = supp_table$pathway,
  `Mean (WT)`                 = sig3(supp_table$mean_WT),
  `Mean (Cd14-/-)`            = sig3(supp_table$mean_KO),
  `Difference (Cd14-/- - WT)` = sig3(supp_table$obs_diff_KO_vs_WT),
  `Exact permutation p`       = round(supp_table$p_value_exact_perm, 3),
  `BH-adjusted p`             = round(supp_table$p_adj_BH, 3),
  `Direction concordance (%)` = round(supp_table$direction_concordance * 100, 1),
  check.names = FALSE
)

write.csv(supp_table_out, "results/SupplementaryTable_pathway_permutationTest.csv", row.names = FALSE)

supp_table_out

## NOTE for the supplementary table legend: the BH-adjusted p-values were
## computed across ALL pathways tested, not only the subset shown in the table.
## This must be stated in the legend, otherwise a reader could reasonably (but
## incorrectly) assume BH was applied only within the displayed rows.

## ---------------------------------------------------------------------------
## Part 6. Supplementary figure: per-pathway dot + bar plot
## ---------------------------------------------------------------------------
# Group mean +/- SEM bar with individual mouse points overlaid. This is the plot
# that directly answers the reviewer -- it shows n=3 per genotype instead of one
# pooled bar. Same pathway set and order as the supplementary table above.
#
# Can be run standalone in a fresh session by reloading the saved CSVs:
#   strength_mat       <- as.matrix(read.csv("results/pathway_strength_perSample.csv", row.names = 1))
#   meta_sample        <- read.csv("results/sample_group_map.csv")
#   pathway_results_df <- read.csv("results/pathway_permutationTest_results.csv")

plot_pathway_dotbar <- function(pathway_name) {
  x  <- strength_mat[pathway_name, ]
  df <- data.frame(sample = names(x), value = as.numeric(x), group = group_vec[names(x)])
  df$group <- factor(df$group, levels = c("WT", "KO"))

  summ <- aggregate(value ~ group, df, function(v) c(mean = mean(v), sem = sd(v)/sqrt(length(v))))
  summ <- do.call(data.frame, summ); names(summ) <- c("group", "mean", "sem")

  row <- pathway_results_df[pathway_results_df$pathway == pathway_name, ]
  lab <- sprintf("p=%.1f, concordance=%.0f%%", row$p_value_exact_perm, 100 * row$direction_concordance)

  ggplot() +
    geom_col(data = summ, aes(group, mean, fill = group), width = 0.6, alpha = 0.5, color = "black") +
    geom_errorbar(data = summ, aes(group, ymin = mean - sem, ymax = mean + sem), width = 0.2) +
    geom_jitter(data = df, aes(group, value), width = 0.08, size = 2, shape = 21, fill = "white") +
    scale_fill_manual(values = c(WT = "#4C72B0", KO = "#DD8452")) +
    labs(title = pathway_name, subtitle = lab, x = NULL, y = "Pathway strength") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", plot.subtitle = element_text(size = 8, color = "grey30"))
}

supp_figure <- wrap_plots(lapply(pathways_for_table, plot_pathway_dotbar), ncol = 4)

ggsave("results/SupplementaryFigure_perSample_dotbar.pdf", supp_figure, width = 10, height = 8)
