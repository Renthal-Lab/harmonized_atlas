library(hash)
library(Seurat)
#library(rtracklayer)
library(tibble)
library(dplyr)
library(gridExtra)
library(grid)
library(ggplot2)
library(hdf5r)
#library(Hmisc)
library(biomaRt)
library(rliger)

#functions
hstomm <- function(seurat_human) {
  ## change human gene into mouse
  hsGene=rownames(seurat_human@assays$RNA@counts)
  print(hsGene)
  
  
  genesV2 <- genesV2[!duplicated(genesV2$external_gene_name) & genesV2$mmusculus_homolog_associated_gene_name != "" ,] #keeping the first of duplicate for simplicity
  int_gene<-intersect(hsGene,genesV2$external_gene_name)
  diff_gene<-setdiff(hsGene,genesV2$external_gene_name) #to keep anything not in mouse ensembl
  print(length(int_gene))
  print(length(diff_gene))
  int_gene<-c(int_gene,diff_gene)
  df <- data.frame (external_gene_name  = diff_gene,
                    mmusculus_homolog_associated_gene_name = diff_gene
  ) #to keep anything not on Ensembl
  genesV2<-rbind(genesV2,df)
  rownames(genesV2)<-genesV2$external_gene_name
  length(int_gene)
  cd_hs<-seurat_human@assays$RNA@counts[int_gene,]#raw counts
  
  rownames(cd_hs)<-genesV2[int_gene,"mmusculus_homolog_associated_gene_name"]
  #cd_human <- cd_human[!duplicated(genesV2$mmusculus_homolog_associated_gene_name),]
  temp_meta <- seurat_human@meta.data
  temp_meta$bc <- rownames(temp_meta)
  temp_meta$nFeature_RNA <- NULL
  temp_meta$nCount_RNA <- NULL
  seurat_human <- CreateSeuratObject(cd_hs)
  seurat_human <- NormalizeData(seurat_human)
  seurat_human@meta.data$bc <- rownames(seurat_human@meta.data)
  seurat_human@meta.data <- left_join(seurat_human@meta.data,temp_meta, by = "bc")
  rownames(seurat_human@meta.data) <- seurat_human@meta.data$bc 
  return(seurat_human)
}
human_mouse_ortho<-read.csv("~/harmonized_DRG_TG_atlasv1.3/orthology_conversion/human_mouse_orthos.txt")
seurat_human<-readRDS("/n/scratch3/users/s/shb953/25SEPT2023/drg_data/Jung_2023//human/Seurat.rds")
DefaultAssay(seurat_human)<-"RNA"
seu_converted<-hstomm(seurat_human)

seu_converted@meta.data<-seurat_human@meta.data
#seu_converted$Publication<-"Jung_2023"
#seu_converted$Platform<-"10x"
#seu_converted$Species<-"Human"
#saveRDS(seu_converted,"/n/scratch3/users/s/shb953/25SEPT2023/drg_data/Jung_2023/human/Seurat_biomart.Rds")
