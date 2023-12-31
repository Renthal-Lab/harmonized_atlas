###Skeleton code for our Seurat integration of Trigeminal ganglia sc/snRNA-seq data
###Author: Renthal Lab
###July 5, 2023


####PACKCAGES#####
library(Seurat)
library(tidyverse)
library(grid)

####GLOBAL VARIABLES#####
pca_num <- 30
res <- 1
## minimum Feature = 400
## mitochondrial gene < 10%

minFeature = 400
maxpercent.mt = 10
nfeature_tg_int <- 2000


####DATA INPUT#####
setwd("~")
seurat_ginty_nature_2020 <- readRDS("Ginty_nature_2020/20210517_Sharma_TG_level1_raw_600_1/Seurat_anchored.Rds")
seurat_jerome_natcomm_2020 <- readRDS("Jerome_natcomm_2020/Seurat.Rds")
seurat_renthal_neuron_2022_human <- readRDS("Renthal_neuron_2022/human/20210724_human_TG_level1_final_updated/Seurat_biomart.Rds")
seurat_renthal_neuron_2022_mouse <- readRDS("Renthal_neuron_2022/mouse/20220102_mouse_TG_merge/Seurat.Rds")
seurat_ryba_elife_2019_10x <- readRDS("Ryba_elife_2019/Seurat_10genomics_annotated.Rds")
seurat_ryba_elife_2019_drop <- readRDS("Ryba_elife_2019/Seurat_drop_annotated.Rds")
seurat_ryba_plos_2017 <- readRDS("Ryba_plos_2017/Seurat_annotated.Rds")

####DATA WRANGLING#####
## ginty_2020
seurat_ginty_nature_2020@meta.data$Age <- "P28–42"
seurat_ginty_nature_2020@meta.data$Strain <- "C57Bl/6J"
seurat_ginty_nature_2020@meta.data$Species <- "Mm"
seurat_ginty_nature_2020@meta.data$platform <- "10xsc"
seurat_ginty_nature_2020@meta.data$publication <- "ginty_nature_2020"

## jermore nature communication
seurat_jerome_natcomm_2020@meta.data$Age <- "6-8w"
seurat_jerome_natcomm_2020@meta.data$Strain <- "Swiss"
seurat_jerome_natcomm_2020@meta.data$Sex <- "F"
seurat_jerome_natcomm_2020@meta.data$Species <- "Mm"
seurat_jerome_natcomm_2020@meta.data$platform <- "10xsc"
seurat_jerome_natcomm_2020@meta.data$publication <- "jerome_natcomm_2020"

## Renthal neuron 2022 human
seurat_renthal_neuron_2022_human@meta.data$Species <- "Hs"
seurat_renthal_neuron_2022_human@meta.data$platform <- "10xsn"
seurat_renthal_neuron_2022_human@meta.data$publication <- "renthal_neuron_2022_human"
seurat_renthal_neuron_2022_human@meta.data$Sex <- gsub("[12]","F",seurat_renthal_neuron_2022_human@meta.data$donor)
seurat_renthal_neuron_2022_human@meta.data$Sex <- gsub("[3]","M",seurat_renthal_neuron_2022_human@meta.data$Sex)

## Renthal neuron 2022 mouse 
seurat_renthal_neuron_2022_mouse@meta.data$Age <- "8-12w"
seurat_renthal_neuron_2022_mouse@meta.data$Strain <- "C57Bl/6J"
seurat_renthal_neuron_2022_mouse@meta.data$Species <- "Mm"
seurat_renthal_neuron_2022_mouse@meta.data <- seurat_renthal_neuron_2022_mouse@meta.data %>% dplyr::rename(Sex = sex)
seurat_renthal_neuron_2022_mouse@meta.data$Sex <- gsub("male","M",seurat_renthal_neuron_2022_mouse@meta.data$Sex)
seurat_renthal_neuron_2022_mouse@meta.data$Sex <- gsub("Male","M",seurat_renthal_neuron_2022_mouse@meta.data$Sex)
seurat_renthal_neuron_2022_mouse@meta.data$Sex <- gsub("Female","F",seurat_renthal_neuron_2022_mouse@meta.data$Sex)

## separate Renthal mouse based on platform
seurat_renthal_neuron_2022_mouse_indrops <- seurat_renthal_neuron_2022_mouse %>% subset(orig.ident == "Reference" | orig.ident == "Anchor_indrops")
seurat_renthal_neuron_2022_mouse_10x <- seurat_renthal_neuron_2022_mouse %>% subset(orig.ident == "Anchor_10X")
seurat_renthal_neuron_2022_mouse_indrops$platform <- "indrops_sn"
seurat_renthal_neuron_2022_mouse_10x$platform <- "10xsn"
seurat_renthal_neuron_2022_mouse_indrops@meta.data$publication <- "renthal_neuron_2022_mouse_indrops"
seurat_renthal_neuron_2022_mouse_10x@meta.data$publication <- "renthal_neuron_2022_mouse_10x"

## Ryba elife 2019 10x
seurat_ryba_elife_2019_10x@meta.data$Age <- "6w_or_older"
seurat_ryba_elife_2019_10x@meta.data$Strain <- "C57BL/6NCrl"
seurat_ryba_elife_2019_10x@meta.data$Species <- "Mm"
seurat_ryba_elife_2019_10x@meta.data$platform <- "10xsn"
seurat_ryba_elife_2019_10x@meta.data$publication <- "ryba_elife_2019_10x"

## Ryba elife 2019 drop
seurat_ryba_elife_2019_drop@meta.data$Age <- "6w_or_older"
seurat_ryba_elife_2019_drop@meta.data$Strain <- "C57BL/6NCrl"
seurat_ryba_elife_2019_drop@meta.data$Species <- "Mm"
seurat_ryba_elife_2019_drop@meta.data$platform <- "drop_sn"
seurat_ryba_elife_2019_drop@meta.data$publication <- "ryba_elife_2019_drop"

# seurat ryba plos 2017
seurat_ryba_plos_2017@meta.data$Age <- "3–5w"
seurat_ryba_plos_2017@meta.data$Sex <- "F/M"
seurat_ryba_plos_2017@meta.data$Strain <- "FVB/N"
seurat_ryba_plos_2017@meta.data$Species <- "Mm"
seurat_ryba_plos_2017@meta.data$platform <- "drop_sc"
seurat_ryba_plos_2017@meta.data$publication <- "ryba_plos_2017"



## subset meta.data
seurat_ginty_nature_2020@meta.data <- seurat_ginty_nature_2020@meta.data %>% 
  dplyr::select(c("Age","Strain","Species","subtype","platform","publication","percent.mt","nCount_RNA","nFeature_RNA"))
seurat_jerome_natcomm_2020@meta.data <- seurat_jerome_natcomm_2020@meta.data %>%
  dplyr::select(c("Age","Sex","Strain","Species","platform","publication","percent.mt","nCount_RNA","nFeature_RNA"))
seurat_renthal_neuron_2022_human@meta.data <- seurat_renthal_neuron_2022_human@meta.data %>%
  dplyr::select(c("Sex","Species","subtype","platform","publication","percent.mt","nCount_RNA","nFeature_RNA"))
seurat_renthal_neuron_2022_mouse_indrops@meta.data <- seurat_renthal_neuron_2022_mouse_indrops@meta.data %>%
  dplyr::select(c("Age","Strain","Sex","Species","subtype","platform","publication","percent.mt","nCount_RNA","nFeature_RNA"))
seurat_renthal_neuron_2022_mouse_10x@meta.data <- seurat_renthal_neuron_2022_mouse_10x@meta.data %>%
  dplyr::select(c("Age","Strain","Sex","Species","subtype","platform","publication","percent.mt","nCount_RNA","nFeature_RNA"))
seurat_ryba_elife_2019_10x@meta.data <- seurat_ryba_elife_2019_10x@meta.data %>%
  dplyr::select(c("Age","Strain","Species","subtype","platform","publication","percent.mt","nCount_RNA","nFeature_RNA"))
seurat_ryba_elife_2019_drop@meta.data <- seurat_ryba_elife_2019_drop@meta.data %>%
  dplyr::select(c("Age","Strain","Species","subtype","platform","publication","percent.mt","nCount_RNA","nFeature_RNA"))
seurat_ryba_plos_2017@meta.data <- seurat_ryba_plos_2017@meta.data %>%
  dplyr::select(c("Age","Strain","Sex","Species","subtype","platform","publication","percent.mt","nCount_RNA","nFeature_RNA"))

## subset data based on QC
seurat_renthal_neuron_2022_mouse_indrops <- subset(seurat_renthal_neuron_2022_mouse_indrops, subset = nFeature_RNA > minFeature)
seurat_renthal_neuron_2022_mouse_10x <- subset(seurat_renthal_neuron_2022_mouse_10x, subset = nFeature_RNA > minFeature)
seurat_ryba_plos_2017 <- subset(seurat_ryba_plos_2017, subset = nFeature_RNA > minFeature)

seurat_jerome_natcomm_2020 <- subset(seurat_jerome_natcomm_2020, percent.mt < maxpercent.mt)
seurat_ryba_plos_2017 <- subset(seurat_ryba_plos_2017, percent.mt < maxpercent.mt)




setwd("../20220601_seurat_obj_update/")
saveRDS(seurat_ginty_nature_2020,"seurat_ginty_nature_2020.Rds")
saveRDS(seurat_jerome_natcomm_2020,"seurat_jerome_natcomm_2020.Rds")
saveRDS(seurat_renthal_neuron_2022_human,"seurat_renthal_neuron_2022_human.Rds")
saveRDS(seurat_renthal_neuron_2022_mouse_indrops,"seurat_renthal_neuron_2022_mouse_indrops.Rds")
saveRDS(seurat_renthal_neuron_2022_mouse_10x,"seurat_renthal_neuron_2022_mouse_10x.Rds")
saveRDS(seurat_ryba_elife_2019_10x,"seurat_ryba_elife_2019_10x.Rds")
saveRDS(seurat_ryba_elife_2019_drop,"seurat_ryba_elife_2019_drop.Rds")
saveRDS(seurat_ryba_plos_2017,"seurat_ryba_plos_2017.Rds")


seurat_tg_list <- list("ginty_nature_2020" = seurat_ginty_nature_2020,
                        "jerome_natcomm_2020" = seurat_jerome_natcomm_2020,
                        "renthal_neuron_2022_human" = seurat_renthal_neuron_2022_human,
                        "renthal_neuron_2022_mouse_indrops" = seurat_renthal_neuron_2022_mouse_indrops,
                       "renthal_neuron_2022_mouse_10x" = seurat_renthal_neuron_2022_mouse_10x,
                        "ryba_elife_2019_10x" = seurat_ryba_elife_2019_10x,
                        "ryba_elife_2019_drop" = seurat_ryba_elife_2019_drop,
                        "ryba_plos_2017" = seurat_ryba_plos_2017)


####INTEGRATION STEPS FOLLOW TUTORIAL: https://satijalab.org/seurat/articles/integration_introduction.html
# normalize and identify variable features for each dataset independently
seurat_tg_list <- lapply(X = seurat_tg_list, FUN <- function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeature_tg_int)
})

# select features that are repeatedly variable across datasets for integration. Provide reference
# for integration anchors to be from only scRNA-seq data
features <- SelectIntegrationFeatures(object.list = seurat_tg_list)

tg.anchors <- FindIntegrationAnchors(object.list = seurat_tg_list, anchor.features = features)

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
ElbowPlot(tg.combined, ndims = pca_num)
tg.combined <- RunUMAP(tg.combined, reduction = "pca", dims = 1:pca_num)
tg.combined <- RunTSNE(tg.combined, dims = 1:pca_num, reduction = "pca",check_duplicates = FALSE)
tg.combined <- FindNeighbors(tg.combined, reduction = "pca", dims = 1:pca_num)
tg.combined <- FindClusters(tg.combined, resolution = res)

tg.combined = AddMetaData(tg.combined,Embeddings(tg.combined[["tsne"]]),colnames(Embeddings(tg.combined[["tsne"]])))
tg.combined = AddMetaData(tg.combined,Embeddings(tg.combined[["umap"]]),colnames(Embeddings(tg.combined[["umap"]])))

DimPlot(tg.combined,reduction = "tsne",label = T)
DimPlot(tg.combined,reduction = "umap",label = T)

setwd("../20220601_TG_integration_30_1/")
pdf("tSNE.pdf")
DimPlot(tg.combined, reduction = "tsne",label = T)
dev.off()

pdf("UMAP.pdf")
DimPlot(tg.combined, reduction = "umap",label = T)
dev.off()

pdf("QC_Vln.pdf", width = 16, height = 9)
VlnPlot(tg.combined, features = c("nFeature_RNA"), group.by = "publication", pt.size = 0)
VlnPlot(tg.combined, features = c("nCount_RNA"), group.by = "publication", pt.size = 0)
VlnPlot(tg.combined, features = c("percent.mt"), group.by = "publication", pt.size = 0)
dev.off()


saveRDS(tg.combined,"Seurat.Rds")
meta.data <- tg.combined@meta.data
write.csv(meta.data,"meta.data.csv")

pdf("UMAP_publication.pdf")
DimPlot(tg.combined,group.by = "publication", label = F)
dev.off()

pdf("UMAP_Species.pdf")
DimPlot(tg.combined,group.by = "Species", label = F)
dev.off()

##MARKER GENES
#tg.combined <- readRDS("Seurat.Rds")
DefaultAssay(tg.combined) <- "RNA"
markergene_neuron <- list("cLTMR1" = c("Fam19a4","Th"),
                          "p_cLTMR2" = c("Fam19a4","Th"),
                          "PEP1" = c("Tac1","Gpx3","Cartpt","Calca"),
                          "PEP2" = c("Tac1","Hpca","Trpm8","Trpv1"),
                          "NP1" = c("Cd55","Mrgprd","Lpar3"),
                          "NP2" = c("Mrgpra3"),
                          "NF1" = c("Nefh","Hapln4","Htr3a","Cplx2"),
                          "NF2" = c("Nefh","Hapln4","Pvalb","Calb1"),
                          "NF3" = c("Nefh","Cadps2","Ntrk2","Ntrk1","Ntrk3"),
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
grid.draw(FeaturePlot(tg.combined, c("Snap25","Rbfox3","Sparc","Nefh","Piezo2","Tac1"),reduction='tsne'))
dev.off()

pdf(paste0("FeaturePlot_umap_RNA/FeaturePlot_tSNE_RNA_general.pdf"),heigh=ceiling(3/2)*4,width=8)
grid.draw(FeaturePlot(tg.combined, c("Snap25","Rbfox3","Sparc","Nefh","Piezo2","Tac1"),reduction='umap'))
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

pdf(paste0("VlnPlot_RNA/VlnPlot_RNA_general.pdf"),height = 9,width = 16)
for (j in c("Snap25","Rbfox3","Sparc","Nefh","Piezo2","Tac1")) {
  grid.draw(VlnPlot(tg.combined,j,group.by = "seurat_clusters",pt.size = 0))
}
dev.off()

###DIFFERENTIAL EXPRESSION
cells.select=tg.combined@meta.data%>%
  rownames_to_column('V1')%>%
  sample_frac(0.15)%>%
  pull(V1)

print(dim(tg.combined))
tg.combined <- subset(tg.combined,cells=cells.select)
print(dim(tg.combined))

cluster_gene <- FindAllMarkers(tg.combined)
write.csv(cluster_gene,"clusterMarkers_RNA_0.15.csv")



