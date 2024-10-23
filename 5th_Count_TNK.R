### 0. set-up
rm(list=ls());gc()
setwd('/media/yschang/T/MY_TNK/')
outdir = './results/Count/'
library(Seurat); library(dplyr); library(ggplot2); library(ggpubr); library(stringr); library(FSA)

# 0.1 get previous work
TNK.sfj <- readRDS('./results/TC_DEG/TNK.sfj')

### 1. Basic metadata
Cases = c("Case_01", "Case_02", "Case_03", "Case_04", "Case_05", 
          "Case_07", "Case_08", "Case_09", "Case_10", "Case_12")
Tissues = c("NL", "Tu")
Smoking_Hx = c("Never", "Cur")
Never_smoker = c("Case_01", "Case_02", "Case_03", "Case_04", "Case_05")
Cur_smoker = c("Case_07", "Case_08", "Case_09", "Case_10", "Case_12")

### 2. Count by patient.id, tissue.id and subcell_modi.id
TNK.sfj <- SetIdent(TNK.sfj, value = "patient.id")
sink(paste0(outdir, 'TC_cnt_by_case.txt'), split = TRUE)
for (idx in levels(TNK.sfj)) {
  a <- subset(TNK.sfj, idents = idx)
  a <- SetIdent(a, value = "tissue.id")
  for (idy in levels(a)) {
    b <- subset(a, idents = idy)
    local_TC.cnt <- ncol(b)
    b <- SetIdent(b, value = "subcell_modi.id")
    for (idz in levels(b)) {
      c <- subset(b, idents = idz)
      d <- ncol(c)/local_TC.cnt*100
      cat(paste0(idx,":",idy,":",idz,":", round(d, 3), "\n"))
    }
  }
}
sink()

### 3. Add metadata
df_01 <- read.table(paste0(outdir,'TC_cnt_by_case.txt'), sep = ":")

colnames(df_01)[1] <- "patient.id"
colnames(df_01)[2] <- 'Tissue.id'
colnames(df_01)[3] <- "subcell_modi.id"
colnames(df_01)[4] <- "Percentile"

df_01 = df_01 %>% mutate(Smoking_Hx = (ifelse(patient.id == "Case_01"|patient.id == "Case_02"|patient.id == "Case_03"|patient.id == "Case_04"|patient.id == "Case_05", "Never", "Cur")))

df_01$Smoking_Hx <- factor(df_01$Smoking_Hx,
                             levels = c("Never", "Cur"),
                             ordered = TRUE)


### 4. Graph and statistics
my_comparison_A = list(c('Never', 'Cur'))
my_comparison_B = list(c('NL', 'Tu'))

# 4.1 BY tissue
for (idx in unique(df_01$subcell_modi.id)) {
  a = subset(df_01, subcell_modi.id == idx)
  ggboxplot(a, x="Tissue.id", y="Percentile", 
            color = "black", fill="Tissue.id",
            xlab = "Group", ylab = 'Percentile', 
            font.label = list(size = 1, face = "plain"), 
            palette = c("#D2D2D2", "#A0A0A0",  '#6E6E6E'), add="jitter") + 
    theme_bw()+
    ylim(0, 35)+
    theme(text = element_text(size = 7.5)) +
    ggtitle(idx)+
    # stat_compare_means(comparisons = my_comparison_A, method = "wilcox.test", paired = TRUE, label.y = 25,  size =2.5) +
    stat_compare_means(comparisons = my_comparison_B, method = "wilcox.test", paired = TRUE, label.y = 30,  size =2.5) 
  # stat_compare_means(method = "kruskal.test", label.y=40)
  ggsave(paste0(outdir,'by_tissue/' , idx, ".jpeg"), width = 5, height = 4, dpi = 300)
}

# 4.2 By smoking
for (idx in unique(df_01$subcell_modi.id)) {
  a = subset(df_01, subcell_modi.id == idx)
  ggboxplot(a, x="Smoking_Hx", y="Percentile", 
            color = "black", fill="Smoking_Hx",
            xlab = "Group", ylab = 'Percentile', 
            font.label = list(size = 1, face = "plain"), 
            palette = c("#D2D2D2", "#A0A0A0",  '#6E6E6E'), add="jitter") + 
    theme_bw()+
    ylim(0, 35)+
    theme(text = element_text(size = 7.5)) +
    ggtitle(idx)+
stat_compare_means(comparisons = my_comparison_A, method = "wilcox.test", paired = TRUE, label.y = 30,  size =2.5) 
  ggsave(paste0(outdir,'by_smoking/' , idx, ".jpeg"), width = 5, height = 4, dpi = 300)
}

# 4.3 According to tissue separated from smoking status.
for (idx in unique(df_01$subcell_modi.id)) {
  a = subset(df_01, subcell_modi.id == idx)
  ggboxplot(a, x="Tissue.id", y="Percentile", 
          color = "black", fill="Tissue.id",
          xlab = "Group", ylab = 'Percentile', 
          facet.by = "Smoking_Hx",
          font.label = list(size = 1, face = "plain"), 
          palette = c("#D2D2D2", "#A0A0A0",  '#6E6E6E'), add="jitter") + 
  theme_bw()+
  ylim(0, 35)+
  theme(text = element_text(size = 7.5)) +
  ggtitle(idx)+
  stat_compare_means(comparisons = my_comparison_B, method = "wilcox.test", paired = TRUE, label.y = 30,  size =2.5) 
# stat_compare_means(method = "kruskal.test", label.y=40)
ggsave(paste0(outdir, 'Faucet_by_smoking/',idx, ".jpeg"), width = 5, height = 4, dpi = 300)
}

# 4.3.1 Post-Hoc
# dunnTest(Percentile ~ Group, data=a, method="bonferroni")

# 4.4 According to smoking separated from tissue status. 
for (idx in unique(df_01$subcell_modi.id)) {
  a = subset(df_01, subcell_modi.id == idx)
  ggboxplot(a, x="Smoking_Hx", y="Percentile", 
            color = "black", fill="Smoking_Hx",
            xlab = "Group", ylab = 'Percentile', 
            facet.by = "Tissue.id",
            font.label = list(size = 1, face = "plain"), 
            palette = c("#D2D2D2", "#A0A0A0",  '#6E6E6E'), add="jitter") + 
    theme_bw()+
    ylim(0, 45)+
    theme(text = element_text(size = 7.5)) +
    ggtitle(idx)+
    stat_compare_means(comparisons = my_comparison_A, method = "wilcox.test", paired = TRUE, label.y = 25,  size =2.5)
ggsave(paste0(outdir, 'Faucet_by_tissue/', idx, ".jpeg"), width = 5, height = 4, dpi = 300)
}
