library(Seurat)
library(TeachingDemos)
library(SeuratWrappers)
library(grid)

wd='/n/scratch3/users/s/shb953/25SEPT2023////tg_subtype_meta/fastmnn_integration/'
dir = paste0(wd,'/',format(Sys.time(),"%Y%m%d"),"_fastMNN_level1_dropseqBD_only_subsetted_400")
dir.create(dir)
setwd(dir)

txtStart("console.txt")
#read in data as Seurat object
seurat_obj<-readRDS("/n/scratch3/users/s/shb953/25SEPT2023/tg_subtype_meta/Seurat_merged.Rds")
DefaultAssay(seurat_obj)<-"RNA"
print("calculate mt %")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^[Mm][Tt]-")
pdf("QC.pdf")
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3,group.by = "dataset")
dev.off()
table(seurat_obj$dataset)
print("subsetting")
seurat_obj <- subset(seurat_obj, subset =  nFeature_RNA>200 & percent.mt < 10)
#seurat_mat
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- RunFastMNN(object.list = SplitObject(seurat_obj, split.by = "dataset"))
seurat_obj <- RunUMAP(seurat_obj, reduction = "mnn", dims = 1:42)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "mnn", dims = 1:42)
seurat_obj <- FindClusters(seurat_obj)
pdf("integration.pdf")
DimPlot(seurat_obj, group.by = c("dataset"),raster = F)
dev.off()
saveRDS(seurat_obj,"fastmnn.Rds")

#seurat_obj <- readRDS("Seurat.Rds")
DefaultAssay(seurat_obj) <- "RNA"

markergene_neuron <- list("cLTMR1" = c("Fam19a4","Th","Gfra2","Pou4f2"),
                          "p_cLTMR2" = c("Fam19a4","Th"),
                          "PEP1" = c("Tac1","Gpx3","Cartpt","Sstr2","Adra2a","Oprk1"),
                          "PEP2" = c("Tac1","Hpca","Trpm8","Trpv1","Rxfp1"),
                          "NP1" = c("Cd55","Mrgprd","Lpar3","Scn11a","Mrgprb4"),
                          "NP2" = c("Mrgpra3"),
                          "NF1" = c("Nefh","Hapln4","Htr3a","Cplx2"),
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
                             "B_cell" = c("Cd79a","Cd79b","Ighd","Ighm","Ms4a1"),
                             "Neutropil" = c("S100a8","S100a9","Retnlg"),
                             "Endothelia" = c("Pecam1","Cldn5","Egfl7"),
                             "Pericyte" = c("Pdgfrb","Notch3","Kcnj8"),
                             "Fibroblast" = c("Dcn","Pdgfra","Mgp"),
                             "T-Cell" = c("Cd3e","Cd8a","Ccr7","Cd4"),
                             "DC" = c("Fcer1a","Cst3"),
                             "NK" = c("Nkg7")
)

markergene <- c(markergene_neuron,markergene_nonneuron)
dir.create('FeaturePlot_tnse_RNA')
dir.create('FeaturePlot_umap_RNA')
index=0


for (i in names(markergene))
{
  index=index+1
  temp_length=sum(markergene[[i]] %in% rownames(seurat_obj))
  if(temp_length==0){
    next()
  }
  print(i)

  pdf(paste0("FeaturePlot_umap_RNA/FeaturePlot_UMAP_RNA_",i,".pdf"),heigh=ceiling(temp_length/2)*4,width=8)
  grid.draw(FeaturePlot(seurat_obj, markergene[[i]],reduction='umap', order = T,raster=F))
  dev.off()
}

pdf(paste0("FeaturePlot_UMAP_RNA_general.pdf"),heigh=ceiling(3/2)*4,width=8)
grid.draw(FeaturePlot(seurat_obj, c("Snap25","Rbfox3","Sparc","Mpz"),reduction='umap'))
dev.off()




txtStop()
