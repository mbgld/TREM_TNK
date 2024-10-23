### 0. set-up
rm(list = ls()); gc()
library(Seurat); 
library(dplyr)
outdir = './results/'; # input.dir = './datasets/'

### 1.0 get previous work
YS_Itg.sfj <- readRDS("/media/yschang/V/Combi/Results/4_Assign/Sim_2_9/YS_Itg.sfj")
MY.sfj <- readRDS('/media/yschang/V/Specific_cluster/1_MY/Results/MY.sfj'); MY.sfj <- SetIdent(MY.sfj, value = "subcell.id")
NK_T.suj <- readRDS('/media/yschang/V/Specific_cluster/4_NKT/Results/NK_T.suj')

### 1.1 Check
DimPlot(YS_Itg.sfj, reduction = "umap", label = T)
DimPlot(MY.sfj, reduction = "umap", label = T)
DimPlot(NK_T.suj, reduction = "umap", label = T)

### 2.0 Combination and re-integration
MY_TNK.combined <- merge(MY.sfj, y = NK_T.suj, add.cell.ids = c("MY", "TNK"))
MY_TNK.list <- SplitObject(MY_TNK.combined, split.by = 'orig.ident')
MY_TNK.list <- lapply(MY_TNK.list, FUN = function(x){
  x<-NormalizeData(x)
  x<-FindVariableFeatures(x, selection.method = 'vst', nfeatures=2000)
})

# 2.1 Find Anchor
MY_TNK.anchors <- FindIntegrationAnchors(object.list = MY_TNK.list, dims=1:30)

# 2.2 Integration
MY_TNK.sbj <- IntegrateData(anchorset=MY_TNK.anchors, dims=1:30)
DefaultAssay(MY_TNK.sbj) <- 'integrated'; rm(YS_dataset.list, MY_TNK.anchors)

### 3. Scaling: .sbj starting point
MY_TNK.ssj <- ScaleData(MY_TNK.sbj, features = rownames(MY_TNK.sbj))

### 4. Linear dimensional Reduction
MY_TNK.spj <- RunPCA(MY_TNK.ssj, npcs = 30, verbose = FALSE); print(MY_TNK.spj[["pca"]], dims = 1:5, nfeatures=5)
# Visualize data
VizDimLoadings(MY_TNK.spj, dims = 1:2, reduction ="pca"); DimPlot(MY_TNK.spj, reduction = 'pca'); DimHeatmap(MY_TNK.spj, dims = 1, cells = 500, balanced =TRUE); DimHeatmap(MY_TNK.spj, dims = 1:15, cells = 500, balanced =TRUE)

### 5.  NON-LINEAR DIMENSION REDUCTION
ElbowPlot(MY_TNK.spj)
MY_TNK.sjj <- JackStraw(MY_TNK.spj, num.replicate = 100)
MY_TNK.sjj <- ScoreJackStraw(MY_TNK.sjj, dims = 1:20)
JackStrawPlot(MY_TNK.sjj, dims = 1:20); rm(MY_TNK.spj)
saveRDS(MY_TNK.sjj, paste0(outdir, 'MY_TNK.sjj'))

### 6. Determine the number of cluster
MY_TNK.sjj <- readRDS(paste0(outdir, 'MY_TNK.sjj'))
MY_TNK.snj <- FindNeighbors(MY_TNK.sjj, dims = 1:30) 
MY_TNK.snj <- FindClusters(MY_TNK.snj, resolution = 2.0)

### 7. Determine distance between clusters using UMAP
MY_TNK.suj <- RunUMAP(MY_TNK.snj, reduction = 'pca',  dims = 1:30)

### 8. Visualize
# plot 1
jpeg(filename = paste0(outdir, "UMAP_Itg_origin.jpeg"), width = 4000, height=3000, res = 600)
DimPlot(MY_TNK.suj, group.by= 'orig.ident', label = T, label.size=3, repel =T)
dev.off()

# plot 2
jpeg(filename = paste0(outdir, "UMAP_Itg_cluster.jpeg"), width = 4000, height=3000, res = 600)
DimPlot(MY_TNK.suj, label = T, repel = T)
dev.off()

# plot 3
jpeg(filename = paste0(outdir, "PCA_Plot_Itg.jpeg"), width = 4000, height=3000, res = 600)
DimPlot(MY_TNK.suj, reduction = 'pca')
dev.off()

# plot 4
MY_TNK.suj <- SetIdent(MY_TNK.suj, value = "subcell.id")
jpeg(filename = paste0(outdir, "UMAP_Itg_label.jpeg"), width = 4000, height=3000, res = 600)
DimPlot(MY_TNK.suj, label = T, repel = T)
dev.off()

### 9. save
saveRDS(MY_TNK.sjj, paste0(outdir, 'MY_TNK.sjj'))
saveRDS(MY_TNK.suj, paste0(outdir, 'MY_TNK.suj'))