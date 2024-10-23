### 0. set-up
rm(list = ls()); gc()
library(Seurat);library(ggplot2)
setwd('/media/yschang/T/MY_TNK/')
outdir = './results/TC_DEG/'

### 1. Prep
# 1.1 get previous work
TNK.suj <- readRDS('./results/TC_DEG/TNK.suj')

### 2. Rename
# 2.1 Set active idents
TNK.suj <- SetIdent(TNK.suj, value = "subcell.id")

# 2.2 Rename clusters
TNK.suj <- RenameIdents(TNK.suj,
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

TNK.suj[["subcell_modi.id"]] <- TNK.suj@active.ident

TNK.sfj <- TNK.suj

DimPlot(TNK.sfj, label = T, repel = T)

### 2. Select genes
DefaultAssay(TNK.sfj) <- 'RNA'
# 2.1 Take out marker genes
TNK_All.mk <- read.table(paste0(outdir, "TNK_All.mk"), header = TRUE)
TNK_Top3.mk <- TNK_All.mk %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
TNK_gene_set_1 <- unique(TNK_Top3.mk$gene)

# 2.2 Custom gene sets_2
TC.mk <- c("CD3D", "CD3E", "CD3G")
abTC.mk <- c("TRAC", "TRBC1", "TRBC2")
gdTC.mk <- c("TRGC1","TRGC2", "TRDC")

naiveTC.mk <- c("CCR7", "SELL"); # homing_to_LN.mk <- c("SELL", "CCR7"), IL7R is relatively low compared to memory TC. So removed

CD4.mk <- c("CD4")
Treg <- c("FOXP3", "ICOS")
# Tfh <- c("CXCR5",  "IL21", "BCL6") # ICOS, PDCD1 is non-specific for Treg and Tfh. 제거 필요. Tfh exist only in the LN.

CD8.mk <-c("CD8A", "CD8B")
activatedTC.mk <-c("IL2RA", "CD40LG") # activated & retained_in_LN.mk <- c("CD69")
controlledTC.mk <- c("CD274", "PDCD1", "CTLA4" )
effectorTC.mk <- c("PRF1", "GZMB", "GZMH", "GZMK")

memoryTC.mk<- c("IL7R") # IL7R is relatively higher to naive TC. Proliferation and maturation.
Lung_retension.mk <- c("ITGAL", "ITGAE", "ITGA4", "CD44", "CD69", "CX3CR1") # retained_in_periphery.mk <- c("CD69"). CD69 is also activated marker. So CD69 was removed.

NK.mk<- c("NCAM1","FCGR3A")

MY_gene_set_1 <- c(abTC.mk, gdTC.mk, naiveTC.mk,
                   CD4.mk, Treg,  
                   CD8.mk, activatedTC.mk, controlledTC.mk, effectorTC.mk,
                   memoryTC.mk, Lung_retension.mk, 
                   NK.mk)

###
TNK.sfj@active.ident <- factor(TNK.sfj@active.ident, 
                               levels = c("naive_T", 
                                          "CD4_eff",
                                          "Treg",
                                          "CD8_eff_2",
                                          "CD8_eff_1",
                                          "CD8_Trm_2",
                                          "CD8_Trm_1",
                                          "γδT_2", 
                                          "γδT_1",
                                          "NK"))


# 4.1 Final plot_1 (by DEG)
jpeg(filename = paste0(outdir, "TNK_DotPlot.jpeg"), width = 8000, height =3200, res = 600)
DotPlot(TNK.sfj, features = TNK_gene_set_1, cols=c('black', 'red'), col.min=-2.5, col.max=2.5, dot.scale = 7.5) &
  theme(axis.text.x = element_text(size=10)) +
  RotatedAxis()
dev.off()

# 4.2 Final plot_1 (by MY geneset)
jpeg(filename = paste0(outdir, "MY_TNK_gene_DotPlot.jpeg"), width = 8000, height =3200, res = 600)
DotPlot(TNK.sfj, features = MY_gene_set_1, cols=c('darkblue', 'red'), col.min=-2.5, col.max=2.5, dot.scale = 7.5) &
  theme(axis.text.x = element_text(size=10)) +
  RotatedAxis()
dev.off()

### 5. Save
saveRDS(TNK.sfj, paste0(outdir, 'TNK.sfj'))
