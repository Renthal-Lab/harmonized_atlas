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
library(dplyr)
#library(Hmisc)

# ========= Functions =========


# ========= Initialization =========

pca_num <- NULL
res <- 1
## minimum Feature = 400
## mitochondrial gene < 10%

minFeature = 500
maxpercent.mt = 10
nfeature_tg_int <- 2000


#markers =list(c('Rbfox3','Sparc'),
#    c('Pax2','Ebf2','Pou6f2','Pitx2'),
#    c('Mme','Fbn2','Gnb4','Pde8a'),
#    c('Aox1','Pdgfd','Glis3','Ahnak2'),
#    c('Atf3','Sox11','Sprr1a','Jun'),
#    c('Jund','Igfbp3','Cdkn1a','Flrt3'),
#    c('Sema6a','Atf5')
#         )

# make dirctory to save data
wd='/n/scratch3/users/s/shb953/25SEPT2023/tg_subtype_meta/seurat_integration/'
dir = paste0(wd,'/',format(Sys.time(),"%Y%m%d"),"_level1_integration_RPCA_dropseqBD_only_subsetted_",minFeature,"_res_",res)
dir.create(dir)
setwd(dir)
txtStart("console.txt")
# ========= Step 1 =========
# load data
print(dir)

#seurat_mat<-readRDS("/n/scratch3/users/s/shb953/25SEPT2023/tg_subtype_meta/seurat_integration/20231009_level1_integration_dropseqBD_only_subsetted_400_res_1/Seurat_merged.Rds")
seurat_mat<-readRDS("/n/scratch3/users/s/shb953/25SEPT2023/tg_subtype_meta/Seurat.Rds")
DefaultAssay(seurat_mat)<-"RNA"
#seurat_mat<-readRDS("/n/scratch3/users/s/shb953/25SEPT2023/tg_subtype_meta/seurat_integration/20231004_level1_integration_actualnocellbend_Nguyen_400_res_1//Seurat_merged.Rds")
#seurat_mat$dataset<-paste(seurat_mat$Publication,seurat_mat$Platform,seurat_mat$Species,sep = "_")
#set.seed(111)
#seurat_mat <- seurat_mat[,sample(nrow(seurat_mat@meta.data), 10000)]
#seurat_mat

#
# ========= Step 2 =========

# calculate mitohondrial genes for each cell
print("calculate mt %")
seurat_mat[["percent.mt"]] <- PercentageFeatureSet(seurat_mat, pattern = "^[Mm][Tt]-")
pdf("QC.pdf")
VlnPlot(seurat_mat, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3,group.by = "dataset")
dev.off()
table(seurat_mat$dataset)
print("subsetting")
#seurat_mat <- subset(seurat_mat, subset =  nFeature_RNA>200 & percent.mt < 10)
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

features <- SelectIntegrationFeatures(object.list = seurat_tg_list)
seurat_tg_list <- lapply(X = seurat_tg_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
# select features that are repeatedly variable across datasets for integration
tg.anchors <- FindIntegrationAnchors(object.list = seurat_tg_list, reduction="rpca",
                                     reference = c(1,2,6,9,10),dims=1:pca_num)

# this command creates an 'integrated' data assay
tg.combined <- IntegrateData(anchorset = tg.anchors,dims = 1:pca_num)

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

pct <- tg.combined[["pca"]]@stdev / sum(tg.combined[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
pca_num<-co1
print(co1)

tg.combined <- RunUMAP(tg.combined, reduction = "pca", dims = 1:pca_num)
tg.combined <- RunTSNE(tg.combined, dims = 1:pca_num, reduction = "pca",check_duplicates = FALSE)
tg.combined <- FindNeighbors(tg.combined, reduction = "pca", dims = 1:pca_num)
tg.combined <- FindClusters(tg.combined, resolution = res)

tg.combined = AddMetaData(tg.combined,Embeddings(tg.combined[["tsne"]]),colnames(Embeddings(tg.combined[["tsne"]])))
tg.combined = AddMetaData(tg.combined,Embeddings(tg.combined[["umap"]]),colnames(Embeddings(tg.combined[["umap"]])))

DimPlot(tg.combined,reduction = "tsne",label = T)
DimPlot(tg.combined,reduction = "umap",label = T)

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


saveRDS(tg.combined,"Seurat.Rds")
meta.data <- tg.combined@meta.data
write.csv(meta.data,"meta.data.csv")

pdf("UMAP_dataset.pdf")
DimPlot(tg.combined,group.by = "dataset", label = F,raster=F)
dev.off()



#tg.combined <- readRDS("Seurat.Rds")
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
  grid.draw(FeaturePlot(tg.combined, markergene[[i]],reduction='umap', order = T))
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

