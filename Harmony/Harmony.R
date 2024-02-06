library(Seurat)
library(TeachingDemos)
library(SeuratWrappers)
library(harmony)
library(grid)

wd='/n/scratch3/users/s/shb953/25SEPT2023////tg_subtype_meta/harmony_integration/'
dir = paste0(wd,'/',format(Sys.time(),"%Y%m%d"),"_Harmony_level1_dropseqBD_only_subsetted_400")
dir.create(dir)
setwd(dir)

txtStart("console.txt")
#read in data as Seurat object
seurat_mat<-readRDS("/n/scratch3/users/s/shb953/25SEPT2023/tg_subtype_meta///Seurat.Rds")

DefaultAssay(seurat_mat)<-"RNA"
print("calculate mt %")
seurat_mat[["percent.mt"]] <- PercentageFeatureSet(seurat_mat, pattern = "^[Mm][Tt]-")
pdf("QC.pdf")
VlnPlot(seurat_mat, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3,group.by = "dataset")
dev.off()
table(seurat_mat$dataset)
print("subsetting")
seurat_mat <- subset(seurat_mat, subset =  nFeature_RNA>200 & percent.mt < 10)

seurat_mat <- NormalizeData(seurat_mat) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)

pct <- seurat_mat[["pca"]]@stdev / sum(seurat_mat[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
pca_num<-co1
print(co1)

seurat_mat <- RunHarmony(seurat_mat, group.by.vars = "dataset",reference_values=c("Aubert 2019_10x_Mouse","Sharma 2020_10x_Mouse",
                                                                                  "Nguyen 2017_Drop-seq_Mouse","Liu 2022_BD Rhapsody_Mouse",
                                                                                  "Liu 2022_BD Rhapsody_Mouse"))
seurat_mat <- RunUMAP(seurat_mat, reduction = "harmony", dims = 1:pca_num)
seurat_mat <- FindNeighbors(seurat_mat, reduction = "harmony", dims = 1:pca_num) %>% FindClusters()

saveRDS(seurat_mat,"harmony.Rds")
pdf("integration.pdf")
DimPlot(seurat_mat, group.by = c("dataset"),raster = F)
dev.off()


#seurat_obj <- readRDS("Seurat.Rds")
DefaultAssay(seurat_mat) <- "RNA"

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
                          "Injured" = c("Atf3", "Jun", "Sox11"),
                          "General" = c("Snap25","Rbfox3","Rbfox1","Snap25","Mpz")
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
  temp_length=sum(markergene[[i]] %in% rownames(seurat_mat))
  if(temp_length==0){
    next()
  }
  print(i)

  pdf(paste0("FeaturePlot_umap_RNA/FeaturePlot_UMAP_RNA_",i,".pdf"),heigh=ceiling(temp_length/2)*4,width=8)
  grid.draw(FeaturePlot(seurat_mat, markergene[[i]],reduction='umap', order = T))
  dev.off()
}





txtStop()
