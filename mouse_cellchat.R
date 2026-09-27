library(Seurat)
library(CellChat)
library(patchwork)

CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB)

obj <- readRDS('data/cd14_for_cellchat.rds') 

###### Group-wise analysis for each group: CD14 knock-out, KO; Wild-type, WT ##
ko <- subset(obj, subset = group == 'KO') 
wt <- subset(obj, subset = group == 'WT')

cc.ko <- createCellChat(object = ko, group.by = "annot", assay = "originalexp")
cc.wt <- createCellChat(object = wt, group.by = "annot", assay = "originalexp")

cc.list <- list(ko = cc.ko, wt = cc.wt)

for (n in names(cc.list)) {
  cc.list[[n]]@DB <- CellChatDB.use
  cc.list[[n]] <- subsetData(cc.list[[n]])
  cc.list[[n]] <- identifyOverExpressedGenes(cc.list[[n]])
  cc.list[[n]] <- identifyOverExpressedInteractions(cc.list[[n]])
  cc.list[[n]] <- computeCommunProb(cc.list[[n]], type = "triMean")
  cc.list[[n]] <- filterCommunication(cc.list[[n]], min.cells = 10)
  cc.list[[n]] <- computeCommunProbPathway(cc.list[[n]])
  cc.list[[n]] <- aggregateNet(cc.list[[n]])
}

df.net.list <- lapply(cc.list, subsetCommunication)

dir.create('results', showWarnings = FALSE)
for (n in names(df.net.list)) {
  write.csv(
    df.net.list[[n]],
    file = file.path('results', paste0('cellchat_interactions_', n, '.csv')),
    row.names = FALSE
  )
}

###### Comparison analysis of KO vs WT ##
cellchat <- mergeCellChat(cc.list, add.names = names(cc.list))

## Figure 2G. Heatmap of interation number and strength difference
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

## Figure 2H. Information flow ##
rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE, return.data = TRUE)

## Figure 2I-J. Hierarchy plot ##
pathways.show <- c('CXCL', 'CCL')

for (p in pathways.show) {
  weight.max <- getMaxWeight(cc.list, slot.name = "netP", attribute = p)
  
  par(mfrow = c(1, length(cc.list)), xpd = TRUE)
  
  for (i in seq_along(cc.list)) {
    netVisual_aggregate(
      cc.list[[i]],
      signaling = p,
      layout = "circle",
      edge.weight.max = weight.max[1],
      edge.width.max = 10,
      signaling.name = paste(p, names(cc.list)[i])
    )
  }
}