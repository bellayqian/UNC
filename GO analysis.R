## GO Enrichment Analysis for Cluster 7,12
library(clusterProfiler)
library(org.Mm.eg.db)

################################################################################
## GO Enrichment Analysis for iCM/Fibs
closet_iCM_BP <- enrichGO(gene        = closest_genes_iCM$gene_name,
                          OrgDb         = org.Mm.eg.db,
                          keyType       = 'SYMBOL',
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05)
barplot(closet_iCM_BP)
write.table(closet_iCM_BP@result, "DA_iCM_GO.txt", sep="\t")

closet_Fibs_BP <- enrichGO(gene         = closest_genes_Fibs$gene_name,
                           OrgDb         = org.Mm.eg.db,
                           keyType       = 'SYMBOL',
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.05)
barplot(closet_Fibs_BP)
write.table(closet_Fibs_BP@result, "DA_Fibs_GO.txt", sep="\t")


################################################################################
# Compare Cluster 7, 12 after trajectory
DefaultAssay(multiomics.subset) <- "RNA"
cluster7.12 <- FindMarkers(multiomics.subset, ident.1 = "7", ident.2 = "12")
cluster7.12 <- filter(cluster7.12, p_val_adj < 0.05)

cluster7.new <- filter(cluster7.12, avg_log2FC > 0) %>% top_n(250, avg_log2FC)
cluster12.new <- filter(cluster7.12, avg_log2FC < 0) %>% top_n(-250, avg_log2FC)

Cluster7_BP <- enrichGO(gene          = rownames(cluster7),
                        OrgDb         = org.Mm.eg.db,
                        keyType       = 'SYMBOL',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05)
Cluster12_BP <- enrichGO(gene         = rownames(cluster12.new),
                         OrgDb         = org.Mm.eg.db,
                         keyType       = 'SYMBOL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05)

p1 <- barplot(Cluster7_BP, showCategory=20) 
p2 <- barplot(Cluster12_BP, showCategory=20) 
p1+p2

################################################################################
# GSEA for DORC
trim_dedup <- read.table(file = "./Results/trim_dedup.txt", sep = "\t")
iCM_Linkage_BP <- enrichGO(gene          = trim_dedup$genes,
                           OrgDb         = org.Mm.eg.db,
                           keyType       = 'SYMBOL',
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.05)
# write.table(iCM_Linkage_BP@result, "./Coembedding/Results/iCM_GO.txt", sep="\t")
iCM_Linkage_BP@result$LOG10 <- -log(iCM_Linkage_BP@result$p.adjust)
iCM_top3BP <- top_n(iCM_Linkage_BP@result, n = 5)

pdf(file = "./Results/Final GO/iCM_DORC_GOplot_BP.pdf", width = 10, height = 10)
ggplot(iCM_top3BP, aes(x=LOG10, y=reorder(Description, LOG10))) +
  geom_bar(stat = "identity", width = 0.5, fill="#F8766D") +
  theme_classic()+
  theme(text = element_text(size = 25)) +
  xlab("") +
  ylab("") +
  xlim(0,20)
dev.off()

trim_dedup.coembed <- read.table(file = "./Coembedding/Results/trim_dedup.txt", sep = "\t")
co.iCM_Linkage_BP <- enrichGO(gene       = trim_dedup.coembed$genes,
                              OrgDb         = org.Mm.eg.db,
                              keyType       = 'SYMBOL',
                              ont           = "BP",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.01,
                              qvalueCutoff  = 0.05)

# write.table(co.iCM_Linkage_BP@result, "./Coembedding/Results/coembed_iCM_GO.txt", sep="\t")

co.iCM_Linkage_BP@result$LOG10 <- -log(co.iCM_Linkage_BP@result$p.adjust)
co_iCM_top3BP <- top_n(co.iCM_Linkage_BP@result, n = 5)

pdf(file = "./Results/Final GO/coiCM_DORC_GOplot_BP.pdf", width = 10, height = 10)
ggplot(co_iCM_top3BP, aes(x=LOG10, y=reorder(Description, LOG10))) +
  geom_bar(stat = "identity", width = 0.5, fill="#F8766D") +
  theme_classic()+
  theme(text = element_text(size = 25)) +
  xlab("") +
  ylab("")
xlim(0,20)
dev.off()

################################################################################
# Compare Cluster 9,10 
DefaultAssay(multiomics) <- "RNA"
cluster9.10 <- FindMarkers(multiomics, ident.1 = "9", ident.2 = "10")
cluster9.10 <- filter(cluster9.10, p_val_adj < 0.05)

cluster9.new <- filter(cluster9.10, avg_log2FC > 0) %>% top_n(250, avg_log2FC)
cluster10.new <- filter(cluster9.10, avg_log2FC < 0) %>% top_n(-250, avg_log2FC)

Cluster9_BP <- enrichGO(gene          = rownames(cluster9.new),
                        OrgDb         = org.Mm.eg.db,
                        keyType       = 'SYMBOL',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05)
Cluster10_BP <- enrichGO(gene         = rownames(cluster10.new),
                         OrgDb         = org.Mm.eg.db,
                         keyType       = 'SYMBOL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05)

Cluster9_BP@result$LOG10 <- -log(Cluster9_BP@result$p.adjust)
Cluster10_BP@result$LOG10 <- -log(Cluster10_BP@result$p.adjust)

cluster9_top3BP <- top_n(Cluster9_BP@result, n = 5)
cluster10_top3BP <- top_n(Cluster10_BP@result, n = 5)

pdf(file = "./Results/Final GO/cluster9GOplot_BP.pdf", width = 10, height = 10)
ggplot(cluster9_top3BP, aes(x=LOG10, y=reorder(Description, LOG10))) +
  geom_bar(stat = "identity", width = 0.5, fill="#F8766D") +
  theme_classic()+
  theme(text = element_text(size = 25)) +
  xlab("") +
  ylab("") + 
  xlim(0,70)
dev.off()

pdf(file = "./Results/Final GO/cluster10GOplot_BP.pdf", width = 10, height = 10)
ggplot(cluster10_top3BP, aes(x=LOG10, y=reorder(Description, LOG10))) + 
  geom_bar(stat = "identity", width = 0.5, fill="#F8766D") +
  theme_classic() +
  theme(text = element_text(size = 25)) +
  xlab("") +
  ylab("") +
  xlim(0,70)
dev.off()


pdf(file = paste0("./Plots/cluster",i,"_9_10GOplot_BP.pdf"), width = 10, height = 10)
p1 <- barplot(Cluster9_BP, showCategory=10) + ggtitle("Cluster 9")
p2 <- barplot(Cluster10_BP, showCategory=10) + ggtitle("Cluster 10")
p1+p2
dev.off()

################################################################################
# enrichGO for each cluster
# i = c(2,12)
i = c(0,3,11)
# i = c(1,4,5,6,7,8)

name <- paste("clusterMarker.", i, sep = "")
markers <- assign(name, dplyr::filter(clusters, clusters$cluster == i))
write.table(markers, file  = "./Results/Final GO/Cluster0311_Marker.txt", sep = "\t")

Group.GO.BP <- enrichGO(gene          = markers$gene,
                        OrgDb         = org.Mm.eg.db,
                        keyType       = 'SYMBOL',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05)
write.table(Group.GO.BP, file  = "./Results/Final GO/Cluster0311_Marker.txt", sep = "\t")

Group.GO.BP@result$LOG10 <- -log(Group.GO.BP@result$p.adjust)
top3BP <- top_n(Group.GO.BP@result, n = 10)

pdf(file = "./Results/Final GO/cluster0311GOplot_BP.pdf", width = 5, height = 5)
ggplot(top3BP, aes(x=LOG10, y=reorder(Description, LOG10))) + 
  geom_bar(stat = "identity", width = 0.7, fill="red3") +
  theme_classic() +
  xlab("-LOG10(Adjusted pvalue)")
dev.off()



#pdf(file = paste0("Cluster",i,"Barplot.pdf"), width = 15, height = 15)
#dev.off()
#gene.Entrez.id <- bitr(cluster$gene,fromType="SYMBOL",toType="ENTREZID",OrgDb='org.Mm.eg.db')
#cluster.name = i
#clutser.profile <- data.frame(gene.Entrez.id, cluster.name)

#gene.Entrez.id <- bitr(all.markers$gene,fromType="SYMBOL",toType="ENTREZID",OrgDb='org.Mm.eg.db')
#gene_char <- as.character(gene.Entrez.id$ENTREZID)

# Barplot for GO_BP analysis
library(tibble)
library(readxl)
library(ggplot2)
library(ggh4x)
library(devtools)
#install_github("kassambara/easyGgplot2")
library(easyGgplot2)
library(tidyverse)

cluster_BP <- read.delim("Clusters.txt",sep = '\t')
cluster_BP <- cluster_BP %>% 
  add_column(Group = if_else(
    .$Description %in% c("actin filament organization",
                         "cardiac muscle tissue development",
                         "muscle system process",
                         "muscle cell differentiation",
                         "muscle tissue development",
                         "muscle organ development",
                         "regulation of actin cytoskeleton organization",
                         "regulation of actin filament-based process",
                         "regulation of muscle system process",
                         "striated muscle tissue development"),"iCM","Others"))
cluster_BP1 <- cluster_BP
cluster_BP1$Group <- replace(cluster_BP$Group, 
                             cluster_BP$Description %in% c("cell-substrate adhesion",
                                                           "extracellular matrix organization",
                                                           "extracellular structure organization",
                                                           "regulation of supramolecular fiber organization"),"Fibroblast")

newBP <- cluster_BP1[, c('Description', 'p.adjust','Group', 'cluster_name','LOG')]
#write.table(newBP, file  = "Final_BP_table.txt")

i = 8
BPtable <- filter(newBP, newBP$cluster_name == i)

p1 <- ggplot(BPtable, aes(LOG, fct_reorder(Description,Group))) + 
  geom_bar(aes(fill=Group), stat='identity') +
  xlab("-LOG ( adjusted p-value)") + ylab("GO Terms") +
  ggtitle(paste0("GO Biological Process plot for Cluster ",i))+
  scale_fill_discrete(breaks = c("iCM", "Others","Fibroblast")) +
  scale_fill_manual(values=c('#00BFC4','#F8766D',"#999999")) +
  theme_classic()+
  theme(axis.text = element_text(size=20, hjust=1)) 
p1

pdf(file = paste0("Cluster",i,"_Barplot_BP.pdf"), width = 12, height = 12)
p1
dev.off()

# Distribution of Group-related GO Terms
pdf(file = "PercentDistribution.pdf", width = 12, height = 12)
ggplot(newBP, aes(cluster_name,fill=Group)) + 
  geom_bar(position="fill") +
  xlab("Cluster Names") + ylab("Percent of Groups") +
  ggtitle("Percent Distribution of GO Terms in clusters")+
  scale_fill_discrete(breaks = c("iCM", "Others","Fibroblast")) +
  scale_fill_manual(values=c('#00BFC4','#F8766D',"#999999")) +
  theme_classic()+
  theme(text = element_text(size=20))
dev.off()
# Expression heatmap for clusters
all.markers %>%
  group_by(cluster) %>%
  top_n(n = -10, wt = (p_val_adj)) -> top10

DefaultAssay(multiomics) <- "SCT"

pdf(file = "Expression_Heatmap_w810.pdf", width = 30, height = 30)
DoHeatmap(multiomics, features = top10$gene)
dev.off()

updated_top10 <- subset(top10, !(cluster %in% c("8","10")))
pdf(file = "Expression_Heatmap_wo810.pdf", width = 30, height = 30)
DoHeatmap(multiomics, features = updated_top10$gene)
dev.off()
