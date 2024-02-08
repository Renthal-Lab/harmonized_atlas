# ========= README =========
# This scipt is used to cluster and visualize scRNA-seq data.
# Step 1: Take the counts table from 10X cellranger output, and generate seurat object
# Step 2: QC, normalize data, find variable genes
# Step 3: PCA, cluster, dimension reduction
# Step 4: generate results including marker gene plots, tSNE, UMAP
# Step 5: save files, find marker genes for each cluster 

library(Seurat)
#library(rtracklayer)
library(tibble)
library(dplyr)
library(gridExtra)
library(grid)
library(ggplot2)
library(hdf5r)
library(TeachingDemos)

# ========= Initialization =========

# Set initial parameters for analysis
res <- 1 # Resolution for clustering
minFeature = 500 # Minimum number of features (genes) to be detected in a cell
maxpercent.mt = 10 # Maximum percentage of mitochondrial genes allowed
nfeature_tg_int <- 2000 # Number of variable features to identify



# Prepare directory for saving data and results
wd='/n/scratch3/users/s/shb953/25SEPT2023/tg_subtype_meta/seurat_integration/'
dir = paste0(wd,'/',format(Sys.time(),"%Y%m%d"),"_level1_integration_",minFeature,"_res_",res)
dir.create(dir) # Create directory
setwd(dir) # Set working directory to newly created directory
txtStart("console.txt") # Start logging console output
# ========= Step 1 =========
# load data
print(dir)
# Load individual datasets as Seurat objects
study1<-readRDS("/n/scratch3/users/s/shb953/25SEPT2023/tg_data/Sharma_2020/Seurat_annotated.Rds")
#repeat for all datasets

# Merge all datasets into a single Seurat object for integrated analysis
seurat_mat<-merge(x=study1,y=c(study2,study3,study4))


#print("Columns of Seurat object")
#colnames(seurat_mat@meta.data)
seurat_mat$dataset<-paste(seurat_mat$Publication,seurat_mat$Platform,seurat_mat$Species,sep = "_")
table(seurat_mat$dataset)

# Save initial merged Seurat object for later reference
saveRDS(seurat_mat,"Seurat_merged.Rds")
#q()
#
# ======== Step 2: Quality Control and Normalization ========
# Calculate the percentage of mitochondrial gene expression per cell as QC metric
print("calculate mt %")
seurat_mat[["percent.mt"]] <- PercentageFeatureSet(seurat_mat, pattern = "^[Mm][Tt]-")

# Plot QC metrics
pdf("QC.pdf")
VlnPlot(seurat_mat, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3,group.by = "dataset")
dev.off()
table(seurat_mat$dataset)

# Filter cells based on QC metrics
print("subsetting")
seurat_mat <- subset(seurat_mat, subset =  nFeature_RNA>minFeature & percent.mt < maxpercent.mt)
#seurat_mat
table(seurat_mat$dataset)
print("splitting the list")

seurat_tg_list <- SplitObject(seurat_mat, split.by = "dataset") 
seurat_tg_list

# normalize and identify variable features for each dataset independently
print("normalizing")
seurat_tg_list <- lapply(X = seurat_tg_list, FUN <- function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeature_tg_int)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = seurat_tg_list)
tg.anchors <- FindIntegrationAnchors(object.list = seurat_tg_list, anchor.features = features,reference = c()) #list datasets for anchoring

# this command creates an 'integrated' data assay
tg.combined <- IntegrateData(anchorset = tg.anchors)

print(paste0("min Feature after QC",min(tg.combined@meta.data$nFeature_RNA)))
print(paste0("min Count after QC",min(tg.combined@meta.data$nCount_RNA)))
print(paste0("max mito after QC",max(tg.combined@meta.data$percent.mt)))

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(tg.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
#tg.combined <- ScaleData(tg.combined, verbose = FALSE)
tg.combined <- ScaleData(tg.combined, features = rownames(tg.combined))
tg.combined <- RunPCA(tg.combined, features = VariableFeatures(tg.combined))

saveRDS(tg.combined,"Seurat.Rds")
pdf("elbow.pdf")
ElbowPlot(tg.combined, ndims = pca_num)
dev.off()

#Picking the number of PCs that cumulatively contribute to ≥ 90% of the standard deviation of gene expression across datasets and
#where the final PC contributes to less than 5% of the standard deviation of gene expression across datasets
pct <- tg.combined[["pca"]]@stdev / sum(tg.combined[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
pca_num<-co1
print(co1)

# Find neighbors and identify clusters
tg.combined <- RunUMAP(tg.combined, reduction = "pca", dims = 1:pca_num)
tg.combined <- RunTSNE(tg.combined, dims = 1:pca_num, reduction = "pca",check_duplicates = FALSE)
tg.combined <- FindNeighbors(tg.combined, reduction = "pca", dims = 1:pca_num)
tg.combined <- FindClusters(tg.combined, resolution = res)

tg.combined = AddMetaData(tg.combined,Embeddings(tg.combined[["tsne"]]),colnames(Embeddings(tg.combined[["tsne"]])))
tg.combined = AddMetaData(tg.combined,Embeddings(tg.combined[["umap"]]),colnames(Embeddings(tg.combined[["umap"]])))

DimPlot(tg.combined,reduction = "tsne",label = T)
DimPlot(tg.combined,reduction = "umap",label = T)

# Generate and save UMAP and t-SNE plots
pdf("tSNE.pdf")
DimPlot(tg.combined, reduction = "tsne",label = T)
dev.off()
pdf("UMAP.pdf")
DimPlot(tg.combined, reduction = "umap",label = T)
dev.off()

pdf("QC_Vln.pdf", width = 16, height = 9)
VlnPlot(tg.combined, features = c("nFeature_RNA"), group.by = "dataset", pt.size = 0)
VlnPlot(tg.combined, features = c("nCount_RNA"), group.by = "dataset", pt.size = 0)
VlnPlot(tg.combined, features = c("percent.mt"), group.by = "dataset", pt.size = 0)
dev.off()

# ======== Step 5: Saving Results and Finding Marker Genes ========
# Save the final Seurat object and metadata
saveRDS(tg.combined,"Seurat.Rds")
meta.data <- tg.combined@meta.data
write.csv(meta.data,"meta.data.csv")

pdf("UMAP_dataset.pdf")
DimPlot(tg.combined,group.by = "dataset", label = F,raster=F)
dev.off()

# Identify and save marker genes for each cluster
DefaultAssay(tg.combined) <- "RNA"
markergene_neuron <- list("cLTMR1" = c("Fam19a4","Th","Gfra2","Pou4f2"),
                          "p_cLTMR2" = c("Fam19a4","Th"),
                          "PEP1" = c("Tac1","Gpx3","Cartpt","Sstr2","Dcn","Adra2a","Oprk1"),
                          "PEP2" = c("Tac1","Hpca","Trpm8","Trpv1","Rxfp1"),
                          "NP1" = c("Cd55","Mrgprd","Lpar3","Scn11a"),
                          "NP2" = c("Mrgpra3","Mrgprb4"),
                          "NF1" = c("Nefh","Hapln4","Htr3a","Cplx2","Bmpr1b","Smr2"),
                          "NF2" = c("Nefh","Hapln4","Pvalb","Calb1"),
                          "NF3" = c("Nefh","Cadps2","Ntrk2"),
                          "SST" = c("Nppb","Sst","Il31ra"),
                          "Abeta_RA" = c("Nefh","Calb1"),
                          "Abeta_field_SA" = c("Nefh","Ntrk3"),
                          "Propcrioceptors" = c("Nefh","Pvalb"),
                          "Adelta" = c("Nefh","Cadps2","Ntrk2"),
                          "Injured" = c("Atf3", "Jun", "Sox11")
)


markergene_nonneuron <- list("Satglia" = c("Apoe","Fabp7","Ednrb"),
                             "schwann_M" = c("Mpz","Mbp"),
                             "schwann_N" = c("Apoe","Scn7a"),
                             "Macrophage" = c("Lyz2","Mrc1","Csf1r","Ccr2"),
                             "B_cell" = c("Cd79a","Cd79b","Ighd","Ighm"),
                             "Neutropil" = c("S100a8","S100a9","Retnlg"),
                             "Endothelia" = c("Pecam1","Cldn5","Egfl7"),
                             "Pericyte" = c("Pdgfrb","Notch3","Kcnj8"),
                             "Fibroblast" = c("Dcn","Pdgfra","Mgp")
)



pdf(paste0("FeaturePlot_UMAP_RNA_general.pdf"),height=ceiling(3/2)*4,width=8)
grid.draw(FeaturePlot(tg.combined, c("Snap25","Rbfox3","Sparc","Mpz","Calca","Retnlg"),reduction='umap',order = T))
dev.off()
#q()

markergene <- c(markergene_neuron,markergene_nonneuron)
dir.create('FeaturePlot_tnse_RNA')
dir.create('FeaturePlot_umap_RNA')
index=0

for (i in names(markergene))
{
  index=index+1
  temp_length=sum(markergene[[i]] %in% rownames(tg.combined))
  if(temp_length==0){
    next()
  }
  print(i)
  pdf(paste0("FeaturePlot_tnse_RNA/FeaturePlot_tSNE_RNA_",i,".pdf"),heigh=ceiling(temp_length/2)*4,width=8)
  grid.draw(FeaturePlot(tg.combined, markergene[[i]],reduction='tsne', order = T))
  dev.off()
  
  pdf(paste0("FeaturePlot_umap_RNA/FeaturePlot_UMAP_RNA_",i,".pdf"),heigh=ceiling(temp_length/2)*4,width=8)
  grid.draw(FeaturePlot(tg.combined, markergene[[i]],reduction='umap', raster=F, order = T))
  dev.off()
}

pdf(paste0("FeaturePlot_tnse_RNA/FeaturePlot_tSNE_RNA_general.pdf"),heigh=ceiling(3/2)*4,width=8)
grid.draw(FeaturePlot(tg.combined, c("Snap25","Rbfox3","Sparc"),reduction='tsne'))
dev.off()

pdf(paste0("FeaturePlot_umap_RNA/FeaturePlot_tSNE_RNA_general.pdf"),heigh=ceiling(3/2)*4,width=8)
grid.draw(FeaturePlot(tg.combined, c("Snap25","Rbfox3","Sparc"),reduction='umap'))
dev.off()

dir.create('VlnPlot_RNA')
index=0
for (i in names(markergene))
{
  index=index+1
  temp_length=sum(markergene[[i]] %in% rownames(tg.combined))
  if(temp_length==0){
    next()
  }
  print(i)
  pdf(paste0("VlnPlot_RNA/VlnPlot_RNA_",i,".pdf"),height = 9,width = 16)
  for (j in markergene[[i]][markergene[[i]] %in% rownames(tg.combined)]) {
    grid.draw(VlnPlot(tg.combined,j,group.by = "seurat_clusters",pt.size = 0))
  }
  dev.off()
}



cells.select=tg.combined@meta.data%>%
  rownames_to_column('V1')%>%
  sample_frac(0.4)%>%
  pull(V1)

print(dim(tg.combined))
tg.combined <- subset(tg.combined,cells=cells.select)
print(dim(tg.combined))

cluster_gene <- FindAllMarkers(tg.combined)
write.csv(cluster_gene,"clusterMarkers_RNA.csv")




txtStop()

q()
