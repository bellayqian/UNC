library(Signac)
library(Seurat)
#library(SeuratWrappers)
#library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
library(scales)
library(viridis)
library(tradeSeq)
library(slingshot)
set.seed(1234)
##read in the seurat object
MA.subset<-readRDS("MA_RNA_ATAC_subset.rds")
##perform slingshot analysis on the wnn.umap reduction and set the start cluster to fibroblast
MA.time <- slingshot(Embeddings(MA.subset, "wnn.umap"), clusterLabels = MA.subset$Bio.ident, 
                     start.clus = "Fib_2", stretch = 0)
##generate color pal for color labeling of the cluster 
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}
#cell_colors <- cell_pal(MGT.subset$, brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(MA.subset$Bio.ident, hue_pal())
#convert the slingshot result into the slingshot object
MA.ssd<-as.SlingshotDataSet(MA.time)
#plot the slingshot result
pdf("MA_pesudotime.pdf")
plot(reducedDims(MA.ssd), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(MA.ssd, lwd = 2, type = "lineages", col = 'black')
dev.off()
#plot each projected trajectories seperately 
nc <- 2
pt <- slingPseudotime(MA.ssd)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100,end = 0.95)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(MA.ssd), col = colors, pch = 16, cex = 0.5, main = i)
  lines(MA.ssd, lwd = 2, col = 'black', type = 'lineages')
}