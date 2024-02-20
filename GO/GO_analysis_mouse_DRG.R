library(topGO)
library(Seurat)
library(tidyverse)
library(org.Mm.eg.db)

setwd("~")

## load in marker gene and seurat object
marker_mouse <- data.table::fread("mouse_markers_DRG.csv")
mouse_seurat <- readRDS("seurat_final_fixed.Rds")

DefaultAssay(mouse_seurat) <- "RNA"

## subset marker genes based on log2FC & p_val_adj
marker_mouse <- subset(marker_mouse, p_val_adj < 0.05 & avg_log2FC > 0.5)

## set background gene to be all the genes in the seurat object
mouse.gene <- rownames(mouse_seurat)

## calculate how many marker genes in each celltype
subtype_m <- unique(marker_mouse$cluster)
nmarker <- data.frame(marker_mouse %>% group_by(cluster) %>% summarise(n = n()))

## take at most 100 marker genes in each celltype
topmarker_m <- NULL
for (i in subtype_m) {
  subtemp <- marker_mouse%>% subset(cluster == i) %>% top_n(100,avg_log2FC)
  topmarker_m <- rbind(topmarker_m,subtemp)
}

GO_result <- NULL

## run GO analysis
for(i in subtype_m){
  submarker_mouse <- topmarker_m %>% subset(cluster == i)
  goi <- submarker_mouse$gene
  geneList <- factor(as.integer(mouse.gene %in% goi))
  names(geneList) <- mouse.gene
  mouseGO <- new("topGOdata", ontology = "BP", allGenes = geneList, 
                 nodeSize = 5, annot = annFUN.org,  mapping = "org.Mm.eg.db",
                 ID = "symbol")
  #
  #Annotated : number of genes in org.Mm.eg.db which are annotated with the GO-term.
  #Significant : number of genes belonging to your input which are annotated with the GO-term.
  #Expected : show an estimate of the number of genes a node of size 
  #           Annotated would have if the significant genes were to be 
  #           randomly selected from the gene universe.
  #
  print(paste("subtype:",i))
  print(paste("Expressing genes:",length(geneList)))
  print(paste("Up genes:",sum(geneList==1)))
  resultFisher <- runTest(mouseGO, algorithm = "weight01", statistic = "fisher")
  genTab=GenTable(mouseGO, weight01Fisher = resultFisher, topNodes = 500, numChar = 1000)%>%
    mutate(expGenes=length(geneList),sigGenes=sum(geneList==1))
  genTab=genTab%>%mutate(Enrichment=Significant/Expected)
  genTab$subtype <- i
  GO_result <- rbind(GO_result,genTab)
}

write.csv(GO_result,"GO_mouse_DRG.csv")



