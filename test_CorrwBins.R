#install.packages("pheatmap")
setwd("~/Desktop/Multiomics/")
library(pheatmap)
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(Seurat)
library(GenomeInfoDb)
library(BiocGenerics)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(dplyr)
library(IRanges)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
#library(SeuratWrappers)
#library(monocle3)
library (BSgenome.Mmusculus.UCSC.mm10)
library(clusterProfiler)
library(org.Mm.eg.db)
set.seed(1234)

#### DO NOT USE DROC GENE LIST ################
# To find DORC gene List for a gene
multiomics.subset <- readRDS("./IntermediateFile/multiomics_subset.rds")
iCM_Linkage <- read.table("./Results/iCM_PeakLinkage_.txt", sep = '\t',header = T)

# # Non-duplicate gene list
# library(dplyr)
# genes <- iCM_Linkage$gene
# dup.gene <- as.data.frame(genes)
# dedup.gene <- dup.gene %>% 
#   group_by(genes) %>% 
#   dplyr::summarize(count = n())
# dedup.gene <- dedup.gene[order(dedup.gene$count),]
# rankNum <- rank(-dedup.gene$count, ties.method = "first")
# dedup.gene <- data.frame(dedup.gene,rankNum)
# genelist.order <- dedup.gene[order(dedup.gene$rankNum),]
# DORC_geneList <- genelist.order$genes

###### Find gene expression for DORC gene list ########
# mulSCT <- SCTransform(multiomics.subset)
mulSCT <- readRDS("./IntermediateFile/SCT_multiomics_sub.rds")
multiomics.ssd <- readRDS("./IntermediateFile/multiomics.ssd.rds")

# Extract specific lineage expression level
Wgene.exp <- t(GetAssayData(mulSCT, slot = "data")[,])
# Do lineage 3, so 3rd column
Wgene.exp_df <- cbind(slingPseudotime(multiomics.ssd)[,3], Wgene.exp)
colnames(Wgene.exp_df)[1] <- "pseudotime"
Wgene.exp_df <- as.data.frame(Wgene.exp_df[!is.na(Wgene.exp_df[,1]),])

## This is the cell x psuedo-time gene expression matrix
L3.gene.exp <- Wgene.exp_df
saveRDS(L3.gene.exp, "./IntermediateFile/lineage3.exp.NoDORC.rds")
L3.gene.exp <- readRDS("./IntermediateFile/lineage3.exp.NoDORC.rds")
# Try not to separate them into bins

# saveRDS(L3.gene.exp, "./IntermediateFile/lineage3.GeneExpression.rds")
# L3.gene.exp <- readRDS("./IntermediateFile/lineage3.GeneExpression.rds")

# Use that to extract wanted cell from ChromVar result

# Split according to pseudotime. c(16.5,12,8,4,0)
# bin <- cut(L3.gene.exp$pseudotime, breaks=c(16.5,8,4,2,0),
#            labels=c("Bin1","Bin2","Bin3","Bin4"),
#            right=FALSE)
# L3.gene.exp2 <- data.frame(bin,L3.gene.exp)
L3.bin1 <- as.data.frame(t(L3.gene.exp))
# L3.bin1 <- rownames(L3.gene.exp2[L3.gene.exp2$bin=="Bin1",])
# L3.bin2 <- rownames(L3.gene.exp2[L3.gene.exp2$bin=="Bin2",])
# L3.bin3 <- rownames(L3.gene.exp2[L3.gene.exp2$bin=="Bin3",])
# L3.bin4 <- rownames(L3.gene.exp2[L3.gene.exp2$bin=="Bin4",])

## Heatmap showing Differential Expression

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(x = JASPAR2020,
                    opts = list(collection = "CORE", 
                                tax_group = 'vertebrates', all_versions = FALSE))

# Add motif information
multiomics.subset <- AddMotifs(object = multiomics.subset,
                               genome = BSgenome.Mmusculus.UCSC.mm10,
                               pfm = pfm)
saveRDS(multiomics.subset, "./IntermediateFile/multiomics_AddMotif")
mul.addMotif <- readRDS("./IntermediateFile/multiomics_AddMotif")

## Finding overrepresented motifs
da_peaks <- FindMarkers(object = mul.addMotif,
                        ident.1 = c("12","2"),ident.2 = c("3","11"),
                        latent.vars = 'nCount_peaks')

# # Get top differentially accessible peaks
# top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.05, ])

# Test enrichment
# enriched.motifs <- FindMotifs(object = mul.addMotif,
#                               features = top.da.peak)
enriched.motifs <- FindMotifs(object = mul.addMotif,
                              features = rownames(da_peaks)) 
# This contains all the motif and its names
saveRDS(enriched.motifs, "./IntermediateFile/enrichedMotifs.rds")
enriched.motifs <- readRDS("./IntermediateFile/enrichedMotifs.rds")

# Computing motif activities
motif.activity <- RunChromVAR(
  object = mul.addMotif,
  genome = BSgenome.Mmusculus.UCSC.mm10)

DefaultAssay(motif.activity) <- 'chromvar'
Motifs <- motif.activity[['chromvar']]@data #Motif Enrichment score for all cells

saveRDS(Motifs, "./IntermediateFile/motifActivity.rds")
Motifs <- readRDS("./IntermediateFile/motifActivity.rds")

# Find Motif Enrichment Score for each bin
bin1.motif <- as.data.frame(Motifs[, colnames(Motifs) %in% colnames(L3.bin1)])
# bin2.motif <- Motifs[, colnames(Motifs) %in% L3.bin2]
# bin3.motif <- Motifs[, colnames(Motifs) %in% L3.bin3]
# bin4.motif <- Motifs[, colnames(Motifs) %in% L3.bin4]

#findCORRScore <- function (bin1.motif,L3.gene.exp2,biNum){

# bin1.exp <- L3.gene.exp2[L3.gene.exp2$bin==biNum,] # With bin
# bin1.exp <- subset(bin1.exp, select=-c(bin,pseudotime))

bin1.exp <- L3.gene.exp
bin1.exp <- as.data.frame(t(bin1.exp))

motif.name <- enriched.motifs[,c("motif","motif.name")]

# test <- bin1.motif[match(row.names(bin1.exp), row.names(bin1.motif)),"motif"]
bin1.exp <- bin1.exp[toupper(rownames(bin1.exp)) 
                     %in% toupper(enriched.motifs$motif.name),]
enrichREF <- enriched.motifs[toupper(enriched.motifs$motif.name) 
                             %in% toupper(rownames(bin1.exp)),]
bin1.motif <- bin1.motif[rownames(bin1.motif) %in% enrichREF$motif,]

#####################################################################
# Find Correlation 
# X-axis: Pearson correlation btw motif enrichment score and expression of its TF.
# Y-axis: log-transform of the motif enrichment p-value calculated using chromVar.
exp <- as.data.frame(bin1.exp)
score <- as.data.frame(bin1.motif)

exp$rowname <- toupper(row.names(exp))
score$rowname <- toupper(row.names(score))

exp.order <- exp[order(exp$rowname),]
enrich.order <- enrichREF[order(enrichREF$motif.name),]
exp.order$motif.name <- enrich.order$motif

row.names(exp.order) <- exp.order$motif.name
exp.order <- exp.order[order(exp.order$motif.name),]
exp.order <- exp.order[,order(colnames(exp.order))]
exp.order <- subset(exp.order, select=-c(motif.name,rowname))

score <- score[order(score$rowname),]
score.order <- score[,order(colnames(score))]
score.order <- subset(score.order, select=-c(rowname))

score.order <- as.data.frame(t(score.order))
exp.order <- as.data.frame(t(exp.order))
saveRDS(score.order, "./IntermediateFile/ScoreMatrix_Whole_preCorr.rds")
saveRDS(exp.order, "./IntermediateFile/ExpMatrix_Whole_preCorr.rds")

cor.score <- cor(score.order,exp.order)

x <- data.frame(Score=diag(cor.score),row.names=rownames(cor.score))
x$motif.name <- enrich.order$motif.name

y <- data.frame(p.val=-log10(enrich.order$p.adjust),row.names=rownames(cor.score))
y$motif.name <- enrich.order$motif.name
xy <- merge(x,y,by.x="motif.name",by.y="motif.name")


# bin1 <- findCORRScore(bin1.motif,L3.bin1,"Bin1")
saveRDS(xy, "./IntermediateFile/prePlot.allBins.rds")
# bin2 <- findCORRScore(bin2.motif,"Bin2")
# bin3 <- findCORRScore(bin3.motif,"Bin3")
# bin4 <- findCORRScore(bin4.motif,"Bin4")
# 
# saveRDS(bin1, "./IntermediateFile/prePlot.bin1.rds")
# saveRDS(bin2, "./IntermediateFile/prePlot.bin2.rds")
# saveRDS(bin3, "./IntermediateFile/prePlot.bin3.rds")
# saveRDS(bin4, "./IntermediateFile/prePlot.bin4.rds")

xy <- readRDS("./IntermediateFile/prePlot.bin1.rds")
xy <- xy[order(-xy$Score),]

xy$label <- c(head(xy$motif.name,10),rep(NA,length.out=413))

# Scatterplot for Correlation
library(ggplot2)
library(ggrepel)
pdf(file = "./Plots/MotifScoreCorrelation.pdf", width = 5, height = 5)
p <- ggplot(xy, aes(x=p.val, y=Score, label = label)) + 
  geom_point(size = 0.8) +
  theme_classic() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey") +
  labs(x = "-LOG10(p.value)", y = "Pearson Correlation Coefficiency") +
  theme(axis.line = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p + geom_text_repel(nudge_x = 0, nudge_y = 0, size = 2,color = "black")
dev.off()

# geom_hline(yintercept = 10, linetype = "dashed", color = "grey") +
#   geom_vline(xintercept = -54, linetype = "dashed", color = "grey") +
#   scale_color_manual(values = setNames(c('red','black'),c(T, F))) +
# scale_fill_gradientn(colours=z2,na.value = "transparent",breaks=c(-0.5,0,0.5)) +


