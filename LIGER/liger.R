##Libraries
library(dplyr)
library(grid)
library(tibble)
library(data.table)
library(Matrix)
library(Seurat)
library(ggplot2)
library(reshape2)
library(cowplot)
library(rliger)

library(TeachingDemos)

dims = 42

# Define the working directory path for saving outputs
wd='/n/scratch3/users/s/shb953/25SEPT2023/tg_subtype_meta/liger_integration/'
dir = paste0(wd,'/',format(Sys.time(),"%Y%m%d"),"k_",dims,"_liger_integrated_louvain_1_select_genes_2000_k42_ref_based")
dir.create(dir)
setwd(dir)
txtStart("console.txt")

# Load a Seurat object from a specified path. This is the level1 object from Seurat CCA
seurat_int<-readRDS("/n/scratch3/users/s/shb953/25SEPT2023/tg_subtype_meta/seurat.Rds")
# Set the default assay to "RNA" for downstream analysis
DefaultAssay(seurat_int)<-"RNA"
print("Converting Seurat objects to LIGER objects")

# Convert the Seurat object to a LIGER object for integrative analysis
seurat_mat = seuratToLiger(seurat_int, combined.seurat = T, assays.use = 'RNA',meta.var="dataset")

print("Normalize")
# Normalize the LIGER object data
seurat_mat <- normalize(seurat_mat)
print("variable features")

# Select genes for integration based on variance, specifying the number of genes and datasets to use
seurat_mat <- selectGenes(seurat_mat,num.genes=2000,datasets.use=c()) #list datasets for anchors


print("scale")
# Scale the data without centering
seurat_mat <- scaleNotCenter(seurat_mat)

# Optimize the alignment using Alternating Least Squares (ALS) with specified dimensions
print("run ALS")
seurat_mat <- optimizeALS(seurat_mat, k = 42)

# Perform quantile normalization across datasets
print("quantile norm")
seurat_mat <- quantile_norm(seurat_mat)

print("find neighbors")
#seurat_mat <- FindNeighbors(seurat_mat, reduction = "iNMF", dims = 1:20)

# Cluster cells using the Louvain algorithm with a specified resolution
seurat_mat <- louvainCluster(seurat_mat, resolution = 1)
print("clusters")
#seurat_mat <- FindClusters(seurat_mat, resolution = 0.4)
# Dimensional reduction and plotting
print("reduction and plotting")
# Run UMAP for dimensionality reduction and visualization
#seurat_mat <- RunUMAP(seurat_mat, dims = 1:ncol(seurat_mat[["iNMF"]]), reduction = "iNMF")
seurat_mat <- runUMAP(seurat_mat, distance = 'cosine', n_neighbors = 30, min_dist = 0.3)
#d1<-DimPlot(seurat_mat, group.by = c("proj", "ident", "subtype"), ncol = 3)
d1<-plotByDatasetAndCluster(seurat_mat, axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = T)
pdf("liger_umap.pdf",height=8,width=12)
#d1
d1[[1]]+d1[[2]]
dev.off()

# Save the processed LIGER object for future use
saveRDS(seurat_mat, "liger.rds")

# Convert the LIGER object back to a Seurat object and save it
drg_liger_seurat<-ligerToSeurat(seurat_mat)
saveRDS(drg_liger_seurat,"liger_seurat.rds")

pdf("liger_markers.pdf",height=8,width=12)

plotGene(seurat_mat, "Rbfox1",axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE,plot.by = "none")
plotGene(seurat_mat, "Snap25",axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE,plot.by = "none")
plotGene(seurat_mat, "Sparc",axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE,plot.by = "none")
plotGene(seurat_mat, "Mpz",axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE,plot.by = "none")
plotGene(seurat_mat, "Tac1",axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE,plot.by = "none")
plotGene(seurat_mat, "Nefh",axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE,plot.by = "none")
plotGene(seurat_mat, "Retnlg",axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE,plot.by = "none")
dev.off()

Rbfox3 <- plotGene(seurat_mat, "Rbfox3",axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE,plot.by = "none")
Snap25 <- plotGene(seurat_mat, "Snap25",axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE,plot.by = "none",scale.by = "none")
Sparc <- plotGene(seurat_mat, "Sparc",axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE,plot.by = "none")

Tac1 <- plotGene(seurat_mat, "Tac1",axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE)
Cd55 <- plotGene(seurat_mat, "Cd55", axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE)
Th <- plotGene(seurat_mat, "Th", axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE)
Atf3 <- plotGene(seurat_mat, "Atf3", axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE)
Nefh <- plotGene(seurat_mat, "Nefh", axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE)
Htr3a <-plotGene(seurat_mat, "Htr3a", axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE)
Cadps2 <- plotGene(seurat_mat, "Cadps2", axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE)
Gpx3 <- plotGene(seurat_mat, "Gpx3", axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE)



d2<-plot_grid(Rbfox3[[1]],Rbfox3[[2]],Rbfox3[[3]],
Tac1[[1]],Tac1[[2]],Tac1[[3]],
Cd55[[1]],Cd55[[2]],Cd55[[3]],
Th[[1]],Th[[2]],Th[[3]],
Atf3[[1]],Atf3[[2]],Atf3[[3]],
Nefh[[1]],Nefh[[2]],Nefh[[3]],
Htr3a[[1]],Htr3a[[2]],Htr3a[[3]],
Cadps2[[1]],Cadps2[[2]],Cadps2[[3]],
Gpx3[[1]],Gpx3[[2]],Gpx3[[3]],
ncol=3)

pdf("marker_genes.pdf",height=36,width=36)
d2
dev.off()

seurat_mat
str(seurat_mat)
#perform differential expression analysis
cluster_runWilcoxon<-runWilcoxon(seurat_mat,compare.method = "clusters")
datasets_runWilcoxon<-runWilcoxon(seurat_mat,compare.method = "datasets")
write.csv(cluster_runWilcoxon,"cluster_runWilcoxon.csv")
write.csv(datasets_runWilcoxon,"datasets_runWilcoxon.csv")

txtStop()
