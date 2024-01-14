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
library(ggplot2)
library(reshape2)
set.seed(1234)

########################## Separate 4 Bins ####################################
## This is the cell x psuedo-time gene expression matrix for Lineage 3
L3.gene.exp <- readRDS("./IntermediateFile/lineage3.exp.NoDORC.rds")

# Split according to pseudotime. c(16.5,12,8,4,0)
bin <- cut(L3.gene.exp$pseudotime, breaks=c(16.5,8,4,2,0),
           labels=c("Bin1","Bin2","Bin3","Bin4"),
           right=FALSE)
L3.gene.exp2 <- data.frame(bin,L3.gene.exp)

L3.bin1 <- rownames(L3.gene.exp2[L3.gene.exp2$bin=="Bin1",])
L3.bin2 <- rownames(L3.gene.exp2[L3.gene.exp2$bin=="Bin2",])
L3.bin3 <- rownames(L3.gene.exp2[L3.gene.exp2$bin=="Bin3",])
L3.bin4 <- rownames(L3.gene.exp2[L3.gene.exp2$bin=="Bin4",])

## Finding overrepresented motifs
# Test enrichment
enriched.motifs <- readRDS("./IntermediateFile/enrichedMotifs.rds")
# Computing motif activities
Motifs <- readRDS("./IntermediateFile/motifActivity.rds")

# Find Motif Enrichment Score for each bin
bin1.motif <- as.data.frame(Motifs[, colnames(Motifs) %in% L3.bin1])
bin2.motif <- as.data.frame(Motifs[, colnames(Motifs) %in% L3.bin2])
bin3.motif <- as.data.frame(Motifs[, colnames(Motifs) %in% L3.bin3])
bin4.motif <- as.data.frame(Motifs[, colnames(Motifs) %in% L3.bin4])

bin1.motif <- bin4.motif

bin1.exp <- L3.gene.exp2[L3.gene.exp2$bin=="Bin4",] # With bin
bin1.exp <- subset(bin1.exp, select=-c(bin,pseudotime))
bin1.exp <- as.data.frame(t(bin1.exp)) # gene x cell matrix

motif.name <- enriched.motifs[,c("motif","motif.name")]

bin1.exp <- bin1.exp[toupper(rownames(bin1.exp)) 
                     %in% toupper(enriched.motifs$motif.name),]
enrichREF <- enriched.motifs[toupper(enriched.motifs$motif.name) 
                             %in% toupper(rownames(bin1.exp)),]
bin1.motif <- bin1.motif[rownames(bin1.motif) %in% enrichREF$motif,]

exp <- as.data.frame(bin1.exp)
score <- as.data.frame(bin1.motif)

exp$rowname <- toupper(row.names(exp))
score$rowname <- toupper(row.names(score))

score <- score[order(score$rowname),]
enrich.order <- enrichREF[order(enrichREF$motif),]
score$motif.name <- toupper(enrich.order$motif.name)

exp.order <- exp[order(exp$rowname),]
row.names(exp.order) <- exp.order$rowname
exp.order <- exp.order[,order(colnames(exp.order))]
exp.order <- subset(exp.order, select=-c(rowname))

row.names(score) <- score$motif.name
score <- score[order(score$motif.name),]
score.order <- score[,order(colnames(score))]
score.order <- subset(score.order, select=-c(motif.name,rowname))

score.order <- as.data.frame(t(score.order))
exp.order <- as.data.frame(t(exp.order))

saveRDS(score.order, paste0("./IntermediateFile/ScoreMatrix_bin4_preCorr.rds"))
saveRDS(exp.order, paste0("./IntermediateFile/ExpMatrix_bin4_preCorr.rds"))

bin1.score <- readRDS("./IntermediateFile/ScoreMatrix_bin1_preCorr.rds")
bin1.exp <- readRDS("./IntermediateFile/ExpMatrix_bin1_preCorr.rds")
bin2.score <- readRDS("./IntermediateFile/ScoreMatrix_bin2_preCorr.rds")
bin2.exp <- readRDS("./IntermediateFile/ExpMatrix_bin2_preCorr.rds")
bin3.score <- readRDS("./IntermediateFile/ScoreMatrix_bin3_preCorr.rds")
bin3.exp <- readRDS("./IntermediateFile/ExpMatrix_bin3_preCorr.rds")
bin4.score <- readRDS("./IntermediateFile/ScoreMatrix_bin4_preCorr.rds")
bin4.exp <- readRDS("./IntermediateFile/ExpMatrix_bin4_preCorr.rds")

# ScoreMatrix <- readRDS("./IntermediateFile/ScoreMatrix_Whole_preCorr.rds")
# ExpMatrix <- readRDS("./IntermediateFile/ExpMatrix_Whole_preCorr.rds")

bin1.score$bin <- "Bin1"
bin1.exp$bin <- "Bin1"
bin2.score$bin <- "Bin2"
bin2.exp$bin <- "Bin2"
bin3.score$bin <- "Bin3"
bin3.exp$bin <- "Bin3"
bin4.score$bin <- "Bin4"
bin4.exp$bin <- "Bin4"

ScoreMatrix <- rbind(bin1.score,bin2.score,bin3.score,bin4.score)
ExpMatrix <- rbind(bin1.exp,bin2.exp,bin3.exp,bin4.exp)

saveRDS(ScoreMatrix,"./IntermediateFile/ScoreMatrix_combineBins.rds")
saveRDS(ExpMatrix,"./IntermediateFile/ExpMatrix_combineBins.rds")

############################### Heatmap #######################################
ScoreMatrix <- readRDS("./IntermediateFile/ScoreMatrix_combineBins.rds")
ExpMatrix <- readRDS("./IntermediateFile/ExpMatrix_combineBins.rds")

binNum <- data.frame(ScoreMatrix$bin)
rownames(binNum) <- rownames(ScoreMatrix)
colnames(binNum) <- c("binNum")
ScoreMatrix <- subset(ScoreMatrix, select=-c(bin))
pdf(file = "./Plots/ScoreMatrix.pdf", width = 5, height = 5)
pheatmap(ScoreMatrix,annotation_row=binNum,cluster_cols = FALSE,cluster_rows = FALSE)
dev.off()


Fibroblast <- c('Postn','Tagln2','Lox','Col3a1','Col8a1','Tcf21',"MS4A6B","COL1A1","FBLN2","ABCA8A")
EC <- c("CDH5","PECAM1","SOX18","FLT1")
Cardiac <- c('Ttn','Dsp','Tnni3','Myl4','Myl7','Ryr2')

df <- as.data.frame(ScoreMatrix[, colnames(Motifs) %in% text])


breaksList = seq(0, 15, by = 1)

# Plots the first heatmap
pheatmap(expressionData[1:10, ], # Plots the first 10 genes of the dataset
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList) # Sets the breaks of the color scale as in breaksList

# Plots the second heatmap with the same color and breaks options
pheatmap(expressionData[20:30, ], # Plots the third 10 genes of the dataset
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList)

hg <- hclust(dist(ScoreMatrix),"complete")
# as.dendrogram(hg) %>% plot(horiz = TRUE)

my_gene_col <- data.frame(cluster = ifelse(test = my_gene_col == 1, yes = "cluster 1", no = "cluster 2"))
my_sample_col <- data.frame(sample = rep(c("tumour", "normal"), c(4,2)))
row.names(my_sample_col) <- colnames(data_subset)

pheatmap(data_subset, annotation_row = my_gene_col, annotation_col = my_sample_col)


melted_cormat <- melt(mydata)
head(melted_cormat)
melted_cormat$new=log2(1+melted_cormat$value)
ggplot(data = melted_cormat, aes(x=variable, y=ID, fill=new)) + 
  geom_tile()

