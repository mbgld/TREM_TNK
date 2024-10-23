### 0. set-up
rm(list = ls()); gc()
setwd('/media/yschang/T1/MY_TNK/')
library(Seurat); library(dplyr)
outdir = './results/TC_DEG/'
TNK.suj <- readRDS('./results/TC_DEG/TNK.suj')

### 1. Define function
MK_merge <- function(x, y){
  z <- merge(x, y, by="row.names", all.x = F, all.y = F)
  z <- z[rev(order(z$avg_log2FC.x, z$power)), ]
  z <- subset(z, select = c('Row.names', 'avg_log2FC.x', 'power', 'pct.1.x', 'pct.2.x', 'p_val', 'p_val_adj'))
}

### 2. Set default assay as RNA for DEG
TNK.suj <- SetIdent(TNK.suj, value = "subcell.id")
DefaultAssay(TNK.suj) <- 'RNA'
TNK.suj$subcell_modi.id <- Idents(TNK.suj)
TNK.suj <- SetIdent(TNK.suj, value = "subcell_modi.id")
TNK.suj <- RenameIdents(TNK.suj,  
                       "γδT_2" = "gd2",
                       "γδT_1" = "gd1")

### 3. FindAllmarkers
TNK_All.mk <- FindAllMarkers(TNK.suj, min.pct = 0.25, only.pos = FALSE, logfc.threshold = 0.25, assay = 'RNA')
write.table(TNK_All.mk, file = paste0(outdir, "/TNK_All_mk.tsv"), quote = F, row.names = FALSE, sep = "\t")

TNK_All.pw <- FindAllMarkers(TNK.suj, min.pct = 0.25, only.pos = FALSE, logfc.threshold = 0.25, assay = 'RNA', test.use = "roc")
TNK_All <- merge(TNK_All.mk, TNK_All.pw, by="row.names", all.x = F, all.y = F)
TNK_All <- subset(TNK_All, select = c('Row.names', 'cluster.x', 'avg_log2FC.x', 'power', 'pct.1.x', 'pct.2.x', 'p_val', 'p_val_adj', 'gene.y'))
TNK_All <- TNK_All[rev(order(TNK_All$cluster.x, TNK_All$avg_log2FC.x, TNK_All$power)), ]
write.table(TNK_All, file = paste0(outdir, "/TNK_All_mp.tsv"), quote = F, sep = "\t")

### 4. Find each cluster markers. 
for (idx in levels(TNK.suj)) {
  a <- FindMarkers(TNK.suj, ident.1 = idx, logfc.threshold = 0.25, assay = 'RNA', min.pct = 0.25, only.pos = FALSE)
  b <- FindMarkers(TNK.suj, ident.1 = idx, logfc.threshold = 0.25, assay = 'RNA', min.pct = 0.25, only.pos = FALSE, test.use = "roc")
  c <- MK_merge(a, b)
  assign(paste0("TNK_",idx,".mp"), c)
  write.table(c, file = paste0(outdir, idx, "_in_TNK_mp.tsv"), quote = F, row.names = FALSE, sep = "\t")
  rm(a, b, c)
}

### 5. Save
saveRDS(TNK.suj, paste0(outdir, 'TNK.suj'))