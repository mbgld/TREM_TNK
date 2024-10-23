###############################################################################
##### Part 0. Clear R environments and get libraries ##########################
###############################################################################
rm(list = ls()); gc();
setwd("/media/yschang/T/MY_TNK/")
library(CellChat); library(patchwork); library(Seurat)
library(future);  future::plan("multisession", workers = 32); options(stringsAsFactors = FALSE);options(future.globals.maxSize = 32 * 1024^3)  # Increase memory limit upto 32 GiB



###############################################################################
##### Part A. Preprocessing dataset and DB ####################################
###############################################################################

### 1. Prep dataset
# 1.1 Get dataset
MY_TNK.suj <- readRDS('./results/MY_TNK.suj')
DefaultAssay(MY_TNK.suj) <- 'RNA'

# 1.2 Rename clusters as 3rd class
MY_TNK.suj <- RenameIdents(MY_TNK.suj,
                        "TCM"="naive_T", 
                        "Activated_CD8_TRM_1" = "CD8_eff_1",
                        "Activated_CD8_TRM_2" = "CD8_eff_2",
                        "CD8_TRM_1" = "CD8_Trm_1",
                        "CD8_TRM_2"= "CD8_Trm_2",
                        "NK_T"="CD4_eff",
                        "Treg"='Treg', 
                        "γδT_1"= "γδT_1", 
                        "γδT_2"="γδT_2",
                        "NK"="NK")

TREM2.cell <- WhichCells(MY_TNK.suj, expression = TREM2 > 2.5)
Idents(object=MY_TNK.suj, cells = TREM2.cell) <- "TREM2(+):Mono-mc"
MY_TNK.suj$subcell.id <- Idents(MY_TNK.suj) # add new subcell Idents (TREM2)

# 1.3 Subset TNK and TREM(+) cells
Group_interest <- c("TREM2(+):Mono-mc", "naive_T", 
                    "CD8_eff_1", "CD8_eff_2",
                    "CD8_Trm_1", "CD8_Trm_2",
                    "CD4_eff", "Treg", 
                    "γδT_1", "γδT_2", "NK")
GI.sfj <- subset(MY_TNK.suj, ident = Group_interest)
table(GI.sfj@active.ident)

### 2. Prep metadata
meta = GI.sfj@meta.data
meta$subcell.id = droplevels(meta$subcell.id, exclude = setdiff(levels(meta$subcell.id),unique(meta$subcell.id)))

### 3. Create Cellchat object
cellchat.obj <- createCellChat(GI.sfj, meta = meta, group.by = 'subcell.id')

cellchat.obj@idents <- factor(cellchat.obj@idents, 
                              levels=c("TREM2(+):Mono-mc",
                                       "naive_T", 
                                       "CD8_eff_1", "CD8_eff_2", "CD8_Trm_1", "CD8_Trm_2",
                                       "CD4_eff", "Treg", "γδT_1", "γδT_2", "NK"))

### 4. Set database
# 4.1 Check CellChatDB
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction) # Show the structure of the database

# 4.2 Select a subset of cellchatDB for cell-cell communication analysis <1. Secreted Signaling, 2. ECM-Receptor, 3. Cell-Cell Contact, 4. Heterodimer, 5. KEGG, 6. Literature, 7. Others 중 택일>
CellChatDB.use <- CellChatDB # 전체 DB 사용
# CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor") # ECM-Receptor DB 사용
outdir = './results/CellChat/ECM/' # DB에 따른 outdir 설정

# 4.3 Set the selected DB in the object
cellchat.obj@DB <- CellChatDB.use

### 5. Preprocessing expression dataset
cellchat.obj <- subsetData(cellchat.obj)
cellchat.obj <- identifyOverExpressedGenes(cellchat.obj)
cellchat.obj <- identifyOverExpressedInteractions(cellchat.obj)

###############################################################################
#### Part B. Inference of cell-cell communications network ####################
###############################################################################

### 1. Compute communication prob & infer cellular communication network
computeAveExpr(cellchat.obj, features = 'TREM2', type='truncatedMean', trim=0.1)
cellchat.obj <- computeCommunProb(cellchat.obj)
cellchat.obj <- filterCommunication(cellchat.obj, min.cells = 10)

### 2. Extract inferred cellular communication network as a df
df.net <- subsetCommunication(cellchat.obj)

### 3. Infer the cell-cell communication at a signaling pathway level
cellchat.obj <- computeCommunProbPathway(cellchat.obj)

### 4. Calculate aggregated cell-cell communication network
## 4.1 Select cell cluster
# 4.1.1 Using all cell cluster
cellchat.obj <- aggregateNet(cellchat.obj)

# 4.1.2 Using TREM2(+):Mono-mc as key cluster
source.cluster = c("TREM2(+):Mono-mc")
target.cluster = c("naive_T", 
                    "CD8_eff_1", "CD8_eff_2", "CD8_Trm_1", "CD8_Trm_2",
                    "CD4_eff", "Treg", "γδT_1", "γδT_2", "NK")
cellchat.obj <- aggregateNet(cellchat.obj,  sources.use = source.cluster, targets.use = target.cluster, remove.isolate = FALSE)

## 4.2 Visualize
groupSize <- as.numeric(table(cellchat.obj@idents))

# 4.2.1 First graph
svg(paste0(outdir,"No_interaction_TREM2_netVisual.svg"), width = 16, height = 10, pointsize = 16)
netVisual_circle(cellchat.obj@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
dev.off()

# 4.2.2 Second graph
svg(paste0(outdir,"Weight_interaction_TREM2_netVisual.svg"), width = 16, height = 10, pointsize = 16)
netVisual_circle(cellchat.obj@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength")
dev.off()

# 4.2.3
mat <- cellchat.obj@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }


###############################################################################
#### Part C. Visualize of cell-cell communication network #####################
###############################################################################
cellchat.obj@netP$pathways
pathways.show <- c("MHC-I") 


### 1. Visualize each signaling pathway
# 1.1 Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat.obj, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# 1.2. Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat.obj, signaling = pathways.show, layout = "circle")

# 1.3. Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat.obj, signaling = pathways.show, layout = "chord")

# 1.4.Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat.obj, signaling = pathways.show, color.heatmap = "Reds")

### 2. Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
netAnalysis_contribution(cellchat.obj, signaling = pathways.show)
pairLR.CXCL <- extractEnrichedLR(cellchat.obj, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair

# 2.1 Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat.obj, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)

# 2.2 Circle plot
netVisual_individual(cellchat.obj, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# 2.3 Chord diagram
netVisual_individual(cellchat.obj, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

### 3. Automatically save the plots of all inferred network
# Access all signaling pathways showing significant communications
pathways.show.all <- cellchat.obj@netP$pathways
# Check the order of cell identity to set suitable vertex.receiver
levels(cellchat.obj@idents)
vertex.receiver = seq(1,4)
# Plot # This code will generate large amout of files!
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
   netVisual(cellchat.obj, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
   gg <- netAnalysis_contribution(cellchat.obj, signaling = pathways.show.all[i])
   ggsave(filename=paste0(outdir, pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

###  4. Visualize cell-cell communication mediated by multiple L-R or signaling pathways
## 4.1 Bubble plot
# 4.1.1 Significant interactions of all signaling 
svg(paste0(outdir,"netVisual_bubble_TREM2.svg"), width = 8, height = 16, pointsize = 12)
netVisual_bubble(cellchat.obj, sources.use = source.cluster, targets.use = target.cluster, remove.isolate = FALSE)
dev.off()
 
# 4.1.2 Significant interactions in certain signaling 
netVisual_bubble(cellchat.obj, sources.use = source.cluster, targets.use = target.cluster, signaling = c("MHC-I"), remove.isolate = FALSE)

# 4.1.3 Significant interactions based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat.obj, signaling = c("MHC-I"))
netVisual_bubble(cellchat.obj, sources.use = source.cluster, targets.use = target.cluster, pairLR.use = pairLR.use, remove.isolate = TRUE)

## 4.2 Chord diagram
# 4.2.1 All significant interactions
# show all the interactions sending from source.cluster
svg(paste0(outdir,"netVisual_chord_gene_from_TREM2.svg"), width = 14, height = 12, pointsize = 12)
netVisual_chord_gene(cellchat.obj, sources.use = source.cluster, targets.use = target.cluster, lab.cex = 0.6, legend.pos.x = 10, legend.pos.y = 50)
dev.off()

# 4.2.2 Show all the interactions received by source.cluster
svg(paste0(outdir,"netVisual_chord_gene_to_TREM2.svg"), width = 14, height = 12, pointsize = 12)
netVisual_chord_gene(cellchat.obj, sources.use = target.cluster, targets.use = source.cluster, lab.cex = 0.6, legend.pos.x = 10, legend.pos.y = 50)
dev.off()

# 4.2.3 Show all the significant interactions associated with certain signaling pathways
netVisual_chord_gene(cellchat.obj, sources.use = source.cluster, targets.use = target.cluster, signaling = c("MHC-I"), legend.pos.x = 8)

# 4.2.4 Show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat.obj, sources.use = source.cluster, targets.use = target.cluster, slot.name = "netP", legend.pos.x = 10)

## 4.3 Violin/dot plot
plotGeneExpression(cellchat.obj, signaling = "MHC-I")
plotGeneExpression(cellchat.obj, signaling = "MHC-I", enriched.only = FALSE)


###############################################################################
#### Part D: Systems analysis of cell-cell communication network ##############
###############################################################################

### 1 Network centrality scores
## 1.1 Compute
cellchat.obj <- netAnalysis_computeCentrality(cellchat.obj, slot.name = "netP")

## 1.2 Visualize
netAnalysis_signalingRole_network(cellchat.obj, signaling = pathways.show, width = 16, height = 4, font.size = 10)

### 2. Signaling role analysis
gg1 <- netAnalysis_signalingRole_scatter(cellchat.obj)
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat.obj, signaling = c("MHC-I"))
gg1 + gg2

### 3. Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.obj, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.obj, pattern = "incoming")
ht1 + ht2

# 3.1 Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat.obj, signaling = c("MHC-I"))

### 4. Identify global communication patterns
## 4.1 Identify and visualize outgoing communication pattern of secreting cells
library(NMF); library(ggalluvial)
selectK(cellchat.obj, pattern = "outgoing")
dev.off()

## 4.2
nPatterns = 3
cellchat.obj <- identifyCommunicationPatterns(cellchat.obj, pattern = "outgoing", k = nPatterns)
dev.off()

# 4.2.1 river plot
netAnalysis_river(cellchat.obj, pattern = "outgoing")

# 4.2.2 dot plot
netAnalysis_dot(cellchat.obj, pattern = "outgoing")

### 5. Identify and visualize incoming communication pattern of target cells
selectK(cellchat.obj, pattern = "incoming")
dev.off()

# 5.1
nPatterns = 4
cellchat.obj <- identifyCommunicationPatterns(cellchat.obj, pattern = "incoming", k = nPatterns)
dev.off()

# 5.2 river plot
netAnalysis_river(cellchat.obj, pattern = "incoming")

# 5.3 dot plot
netAnalysis_dot(cellchat.obj, pattern = "incoming")

### 6. Manifold and classification learning analysis of signaling networks
# 6.1
cellchat.obj <- computeNetSimilarity(cellchat.obj, type = "functional")
cellchat.obj <- netEmbedding(cellchat.obj, type = "functional")
cellchat.obj <- netClustering(cellchat.obj, type = "functional")

# 6.2 Visualization in 2D-space
netVisual_embedding(cellchat.obj, type = "functional", label.size = 3.5)
# netVisual_embeddingZoomIn(cellchat.obj, type = "functional", nCol = 2)

# 6.3 Identify signaling groups based on structure similarity
cellchat.obj <- computeNetSimilarity(cellchat.obj, type = "structural")
cellchat.obj <- netEmbedding(cellchat.obj, type = "structural")
cellchat.obj <- netClustering(cellchat.obj, type = "structural")
# Visualization in 2D-space
netVisual_embedding(cellchat.obj, type = "structural", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat.obj, type = "structural", nCol = 2)
