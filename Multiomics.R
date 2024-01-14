BiocManager::install("BiocManager")
BiocManager::install("EnsDb.Mmusculus.v79")
BiocManager::install("BiocGenerics")
BiocManager::install("GenomeInfoDb")
BiocManager::install("GenomicRanges")
install.packages("Signac")
install.packages("hdf5r")
BiocManager::install("chromVAR")
#BiocManager::install("JASPAR2020")
BiocManager::install("TFBSTools")
#BiocManager::install("motifmatchr")
remotes::install_github('satijalab/seurat-wrappers')
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
BiocManager::install("clusterProfiler")
remove.packages("rlang")
packageVersion("rlang")
install.packages("rlang")
update.packages("rlang")
BiocManager::install("biovizBase")
BiocManager::install("org.Mm.eg.db")
remotes::install_github("satijalab/seurat-wrappers")

library(Signac)
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
#library(monocle3)
library (BSgenome.Mmusculus.UCSC.mm10)
library(clusterProfiler)
library(org.Mm.eg.db)
set.seed(1234)

#setwd("~/R/Multiomics/")
#setwd("/proj/liulab/users/yqian")
setwd("~/Desktop/Multiomics/")

# load the RNA and ATAC data, all in./atac_fragments.tsv.gz
counts <- Read10X_h5("./Inputs/filtered_feature_bc_matrix.h5")
fragpath <- './Inputs/atac_fragments.tsv.gz'

# create a Seurat object containing the RNA data
multiomics <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# get gene annotations for mm10
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "mm10"

saveRDS(annotation, "./IntermediateFile/annotation.rds")
annotation <- readRDS("./IntermediateFile/annotation.rds")

# create ATAC assay and add it to the object
multiomics[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation)
multiomics

saveRDS(multiomics, "./IntermediateFile/multiomics_orig.rds")
multiomics <- readRDS("./IntermediateFile/multiomics_orig.rds")

## Quality Control
# change defaultAssay from RNA to ATAC
DefaultAssay(multiomics) <- "ATAC"

# calculate the strength of the nucleosome signal per cell. 
# computes the ratio of fragments between 147 bp and 294 bp (mononucleosome) to 
# fragments < 147 bp (nucleosome-free)

multiomics <- NucleosomeSignal(multiomics)
# compute TSS enrichment score per cell
multiomics <- TSSEnrichment(multiomics)

# filter out low quality cells
multiomics <- subset(
  x = multiomics,
  subset = nCount_RNA < 30000 &
    nCount_RNA > 1000 &
    nCount_ATAC > 1000 &
    nCount_ATAC < 100000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)

pdf(file = "./Plots/Vlnplot_QC.pdf", width = 15, height = 15)

VlnPlot(
  object = multiomics,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

dev.off()

saveRDS(multiomics, "./IntermediateFile/multiomics_QC.rds")
multiomics <- readRDS("./IntermediateFile/multiomics_QC.rds")

## Gene expression data processing
## Cell-Cycle Scoring and Regression
DefaultAssay(multiomics) <- "RNA"

# normalize the gene expression data using SCTransform, and reduce the dimensionality using PCA.
#multiomics <- SCTransform(multiomics)
# Alternate path
multiomics <- NormalizeData(multiomics)
multiomics <- FindVariableFeatures(multiomics, selection.method = "vst")
multiomics <- ScaleData(multiomics, features = rownames(multiomics))

#Use SCTransform as an alternative to the NormalizeData, FindVariableFeatures, ScaleData workflow. 
####multiomics <- SCTransform(multiomics)

multiomics <- RunPCA(multiomics, features = VariableFeatures(multiomics))

saveRDS(multiomics, "./IntermediateFile/multiomics_CellCycle_RunPCA.rds")
multiomics <- readRDS("./IntermediateFile/multiomics_CellCycle_RunPCA.rds")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

multiomics <- CellCycleScoring(multiomics, s.features = s.genes, 
                               g2m.features = g2m.genes, set.ident = TRUE)

head(multiomics[[]])

# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by phase
multiomics <- RunPCA(multiomics, features = c(s.genes, g2m.genes))
p1 <- DimPlot(multiomics)

multiomics <- ScaleData(multiomics, vars.to.regress = c("S.Score", "G2M.Score"), 
                        features = rownames(multiomics))

saveRDS(multiomics, "./IntermediateFile/multiomics_CellCycle_RregressS_G2M.rds")
multiomics <- readRDS("./IntermediateFile/multiomics_CellCycle_RregressS_G2M.rds")

multiomics <- RunPCA(multiomics)
p2 <- DimPlot(multiomics)

p1+p2

## Joint UMAP visualization
# build a joint neighbor graph using both assays

# RNA analysis
DefaultAssay(multiomics) <- "RNA"
#multiomics <- SCTransform(multiomics, verbose = FALSE) %>% 
multiomics <- RunPCA(multiomics)
multiomics <- RunUMAP(multiomics, dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(multiomics) <- "ATAC"
multiomics <- RunTFIDF(multiomics)
multiomics <- FindTopFeatures(multiomics, min.cutoff = 'q0')
multiomics <- RunSVD(multiomics)
multiomics <- RunUMAP(multiomics, reduction = 'lsi', dims = 2:50, 
                      reduction.name = "umap.atac", reduction.key = "atacUMAP_")

multiomics <- FindMultiModalNeighbors(object = multiomics,reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:50))
multiomics <- RunUMAP(multiomics, nn.name = "weighted.nn", reduction.name = "wnn.umap", 
                      reduction.key = "wnnUMAP_")
multiomics <- FindClusters(multiomics, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

saveRDS(multiomics, "./IntermediateFile/multiomics_wnn_clustered.rds")
multiomics <- readRDS("./IntermediateFile/multiomics_wnn_clustered.rds")

pdf(file = "./Plots/Dimplot_cluster.pdf", width = 12, height = 12)
p1 <- DimPlot(multiomics, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(multiomics, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(multiomics, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

## Create Peaks assay
metadata <- read.csv(file = "./Inputs/per_barcode_metrics.csv", header = TRUE, row.names = 1)
annotation <- readRDS("./IntermediateFile/annotation.rds")

# call peaks using MACS2
peaks <- CallPeaks(multiomics, macs2.path = "./env/bin/macs3")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)

saveRDS(peaks, "./IntermediateFile/peaks_CallPeaks.rds")
peaks <- readRDS("./IntermediateFile/peaks_CallPeaks.rds")

# quantify counts in each peak
macs3_counts <- FeatureMatrix(
  fragments = Fragments(multiomics),
  features = peaks,
  cells = colnames(multiomics))

saveRDS(macs3_counts, "./IntermediateFile/macs3_counts.rds")
macs3_counts <- readRDS("./IntermediateFile/macs3_counts.rds")

# create a new assay using the MACS2 peak set and add it to the Seurat object
multiomics[["peaks"]] <- CreateChromatinAssay(
  counts = macs3_counts,
  fragments = './Inputs/atac_fragments.tsv.gz',
  annotation = annotation)

multiomics

DefaultAssay(multiomics) <- "peaks"

# Add metadata info to multiomics peaks assay
multiomics <- AddMetaData(multiomics, metadata)

saveRDS(multiomics, "./IntermediateFile/multiomics_metadata.rds")
multiomics <- readRDS("./IntermediateFile/multiomics_metadata.rds")

# Make Heatmap for each cluster, specifically showing cluster 9 and 10 are irrelavent
all.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

pdf(file = "./Plots/Heatmap_multiomicsClusters.pdf", width = 12, height = 12)
DoHeatmap(multiomics, features = top10$gene)
dev.off()

# multiomics.subset <- readRDS("./IntermediateFile/multiomics_subset.rds")
# 
# levels(multiomics.subset@active.ident) <- c("Fib_1","iFib_1","iCM_1","Fib_2","iFib_2","iFib_3",
#                                             "iFib_4","iFib_5","iFib_6","Fib_3","iCM_2")
# all.markers
all.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

pdf(file = "./Plots/Heatmap_multiomicsClusters.pdf", width = 12, height = 12)
DoHeatmap(multiomics, features = top10$gene)
dev.off()

multiomics.subset <- subset(x = multiomics, idents = c("0","1","2","3","4","5","6","7","8","11","12"))
saveRDS(multiomics.subset, "./IntermediateFile/multiomics_subset.rds")
multiomics.subset <- readRDS("./IntermediateFile/multiomics_subset.rds")

levels(multiomics.subset@active.ident) <- c("Fib_1","iFib_1","iCM_1","Fib_2","iFib_2","iFib_3",
                               "iFib_4","iFib_5","iFib_6","Fib_3","iCM_2")

pdf(file = "./Plots/Dimplot_subsetCluster.pdf", width = 12, height = 12)
p1 <- DimPlot(multiomics.subset, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(multiomics.subset, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(multiomics.subset, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(file = "./Plots/Dimplot_subsetWNN.pdf", width = 5, height = 5)
DimPlot(multiomics.subset, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
dev.off()

pdf(file = "./Plots/Dimplot_subsetrna.pdf", width = 5, height = 5)
DimPlot(multiomics.subset, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) 
dev.off()

pdf(file = "./Plots/Dimplot_subsetatac.pdf", width = 5, height = 5)
DimPlot(multiomics.subset, reduction = "umap.atac", label.size = 2.5, repel = TRUE)
dev.off()

## Compared between cluster 7 and 12 to determine cluster identity
# DefaultAssay(multiomics.subset) <- "RNA"
# cluster7.12 <- FindMarkers(multiomics.subset, ident.1 = "7", ident.2 = "12")
# cluster7.12 <- filter(cluster7.12, p_val_adj < 0.05)
# cluster7 <- filter(cluster7.12, avg_log2FC > 0)
# cluster12 <- filter(cluster7.12, avg_log2FC < 0)
# Input to GO_Enrichment_Analysis.R

## Find Marker Genes for multiomics dataset
DefaultAssay(multiomics) <- "RNA"
all.markers <- FindAllMarkers(object = multiomics)
saveRDS(all.markers, "./IntermediateFile/all.markers.rds")
all.markers <- readRDS("./IntermediateFile/all.markers.rds")

for (i in 0:12){
  name <- paste("clusterMarker.", i, sep = "")
  markers <- assign(name, dplyr::filter(all.markers, all.markers$cluster ==i))
  write.table(markers, file  = paste0("./Results/Cluster",i,"_Marker.txt"), sep = "\t")
}
# Input to GO_Enrichment_Analysis.R to perform GO analysis

## Find Marker Genes for multiomics.subset dataset
DefaultAssay(multiomics.subset) <- "RNA"
clusters <- FindAllMarkers(multiomics.subset)
# clusters <- filter(clusters, p_val_adj < 0.05)
saveRDS(clusters, "./IntermediateFile/clustersMarkers.rds")
clusters <- readRDS("./IntermediateFile/clustersMarkers.rds")

for (i in 0:12){
  name <- paste("clusterMarker.", i, sep = "")
  markers <- assign(name, dplyr::filter(clusters, clusters$cluster ==i))
  write.table(markers, file  = paste0("./Results/Final GO/Cluster",i,"_Marker.txt"), sep = "\t")
}
# Input to GO_Enrichment_Analysis.R to perform GO analysis


# ###### Does not work ###########
# ## Find differentially accessible peaks between clusters
# DefaultAssay(multiomics.subset) <- 'peaks'
# da_peaks <- FindMarkers(object = multiomics.subset, ident.1 = c("12","2"), ident.2 = c("3"),
#                         min.pct = 0.05, test.use = 'LR', latent.vars = 'atac_peak_region_fragments')
# 
# head(da_peaks)
# 
# plot1 <- VlnPlot(
#   object = multiomics.subset,
#   features = rownames(da_peaks)[1],
#   pt.size = 0.1,
#   idents = c("0","3","11","2","12"))
# 
# plot2 <- FeaturePlot(
#   object = multiomics.subset,
#   features = rownames(da_peaks)[2],
#   pt.size = 0.1)
# 
# plot1 | plot2
# 
# open_iCM <- rownames(da_peaks[da_peaks$avg_log2FC > 0.5, ])
# open_Fibs <- rownames(da_peaks[da_peaks$avg_log2FC < -0.5, ])
# 
# # Empty list for both
# closest_genes_iCM <- ClosestFeature(multiomics.subset, regions = open_iCM)
# closest_genes_Fibs <- ClosestFeature(multiomics.subset, regions = open_Fibs)
# 
# head(closest_genes_iCM)
# head(closest_genes_Fibs)
# ###### Probably does not work ###########

## Find Differentially expressed gene list
DefaultAssay(multiomics.subset) <- 'RNA'
da_motif <- FindMarkers(multiomics.subset,ident.1 = c("2","12"),ident.2 = c("3","11"))

# saveRDS(da_motif, "./IntermediateFile/enriched_Motifs.rds")
# da_motif <- readRDS("./IntermediateFile/enriched_Motifs.rds")

saveRDS(da_motif, "./IntermediateFile/enriched_Motifs.rds")
da_motif <- readRDS("./IntermediateFile/enriched_Motifs.rds")

da_motif <- filter(da_motif, p_val_adj < 0.05)
clusteriCM <- filter(da_motif, avg_log2FC > 0)
clusterFibs <- filter(da_motif, avg_log2FC < 0)

## Find significant differentially accessible gene-peak linkage
## Linking peaks to genes
DefaultAssay(multiomics.subset) <- "peaks"

# first compute the GC content for each peak
multiomics.subset <- RegionStats(multiomics.subset, genome = BSgenome.Mmusculus.UCSC.mm10)

m.Fibs <- LinkPeaks(object = multiomics.subset, peak.assay = "peaks",
                    expression.assay = "RNA",
                    genes.use = rownames(clusterFibs))

m.iCM <- LinkPeaks(object = multiomics.subset, peak.assay = "peaks",
                   expression.assay = "RNA",
                   genes.use = rownames(clusteriCM))

df_Fibs <- as_data_frame(m.Fibs[["peaks"]]@links)
df_iCM <- as_data_frame(m.iCM[["peaks"]]@links)

# write.table(df_Fibs, "./Results/Fibs_PeakLinkage.txt", sep="\t", col.names = NA, row.names = TRUE)
# write.table(df_iCM, "./Results/iCM_PeakLinkage.txt", sep="\t", col.names = NA, row.names = TRUE)

write.table(df_Fibs, "./Results/Fibs_PeakLinkage_.txt", sep="\t", col.names = NA, row.names = TRUE)
write.table(df_iCM, "./Results/iCM_PeakLinkage_.txt", sep="\t", col.names = NA, row.names = TRUE)

# Input to DORC.R to compute DORCs

# Using Joint Analysis significant peak-gene linkage to perform Motif Enrichment Analysis
# iCM.old.Linkage <- read.table("./Results/iCM_PeakLinkage.txt", sep = '\t',header = T)
iCM_Linkage <- read.table("./Results/iCM_PeakLinkage_.txt", sep = '\t',header = T)
# 2022-08-11: Difference between iCM_PeakLinkage_.txt and iCM_PeakLinkage.txt could be FindMarkers for
# _ is between ("2","12") and ("3","11"), while later was between ("2","12") and ("0",3","11"). 

gene <- iCM_Linkage$gene
peakID <- iCM_Linkage$peak
chromosome <- iCM_Linkage$seqnames
start <- iCM_Linkage$start
end <- iCM_Linkage$end
strand <- iCM_Linkage$strand

dup.linkage <- data.frame(gene, peakID, chromosome, start, end, strand)
dedup.linkage <- list()

# trim_dedup is generated from DORC.R line 19
for (i in 1:length(trim_dedup$genes)){
  # print(i)
  for (j in 1: length(dup.linkage$gene)){
    if (trim_dedup[i,1] == dup.linkage[j,1]){
      dedup.linkage <- rbind(dedup.linkage,dup.linkage[j,])
    }
  }
}

write.table(dedup.linkage, "./Coembedding/Results/dedupLinkage.txt", sep="\t", 
            quote=FALSE, col.names = FALSE, row.names = FALSE)

## Output for Homer in Linux
# ssh yqian9@longleaf.unc.edu
# scp dedupLinkage.txt yqian9@longleaf.unc.edu:/nas/longleaf/home/yqian9/R
# 
# cut -f 3,4,5 dedupLinkage.txt > Linkage.bed
# 
# perl /nas/longleaf/home/yqian9/homer/bin/findMotifsGenome.pl /nas/longleaf/home/yqian9/R/Linkage.bed mm10 /nas/longleaf/home/yqian9/R/MotifOutputwparse -size given -preparse
# 
# scp -r yqian9@longleaf.unc.edu:/nas/longleaf/home/yqian9/R/MotifOutputwparse ./
#   
# # Install Homer
#   
# scp configureHomer.pl yqian9@longleaf.unc.edu:/nas/longleaf/home/yqian9/homer
# perl /nas/longleaf/home/yqian9/homer/configureHomer.pl -install
# perl /nas/longleaf/home/yqian9/homer/configureHomer.pl -list
# perl /nas/longleaf/home/yqian9/homer/configureHomer.pl -install mm10
#
# cd /nas/longleaf/home/yqian9/R
# 
# cd /nas/longleaf/home/yqian9/homer/bin

# link peaks to genes
## Ryr2
multiomics.subset <- LinkPeaks(object = multiomics.subset, peak.assay = "peaks",
                               expression.assay = "RNA",
                               genes.use = c("Ryr2"))

pdf(file = "./Plots/CoveragePlot_Ryr2.pdf", width = 12, height = 12)
CoveragePlot(object = multiomics.subset, region = "Ryr2", features = "Ryr2", expression.assay = "RNA",
             extend.upstream = 500, extend.downstream = 10000)
dev.off()

# FeaturePlot for DORCs
DefaultAssay(multiomics.subset) <- "peaks"
multiomics.subset <- RegionStats(multiomics.subset, genome = BSgenome.Mmusculus.UCSC.mm10)
DORCs <- c("Hand2","Atf7","Fos","Foxo1","FOXA1","Sox21","Sox6","Nkx2-5","Smad3","Nkx6-1",
           "Smad4","Mef2c","Tbx5","Gata4")
peakLinkage <- LinkPeaks(object = multiomics.subset, peak.assay = "peaks",
                               expression.assay = "RNA",
                               genes.use = DORCs)

CoveragePlot(object = peakLinkage, region = "Smad3", features = "Smad3", expression.assay = "RNA",
             extend.upstream = 500, extend.downstream = 10000)

pdf(file = "./Plots/FeaturePlot_DORCs.pdf", width = 12, height = 12)
DefaultAssay(multiomics.subset) <- "RNA"
FeaturePlot(multiomics.subset, reduction = "wnn.umap", 
            features = DORCs)
dev.off()


