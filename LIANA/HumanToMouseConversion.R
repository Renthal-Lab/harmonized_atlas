# This code was used to convert any Seurat objects with human gene names in their counts matrix to mouse gene names for the purpose of downstream LR-interaction analysis using LIANA, as mentioned in the Methods of the paper associated with the code

setwd("/n/scratch/users/p/pab1164/pab1164/ShamsLR/TGMeninges_v2")
library(Seurat)
library(tidyverse)
library(tibble)
library(biomaRt)

obj.h <- readRDS("HumanDura.rds")

m.genes=rownames(obj.h@assays$RNA@counts)
human_mart <- biomaRt::useMart(host="https://dec2021.archive.ensembl.org/", "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl") 
mouse_mart <- biomaRt::useMart(host="https://dec2021.archive.ensembl.org/", "ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
human_to_mouse_homologs = biomaRt::getLDS(attributes = c("hgnc_symbol","entrezgene_id","ensembl_gene_id"), filters = "hgnc_symbol", values = m.genes,mart = human_mart,attributesL = c("mgi_symbol","ensembl_gene_id","entrezgene_id"), martL = mouse_mart)  #map human genes to their mouse homologs

#The below section creates a tibble with human genes and then joins it with the human-to-mouse mapping to keep only mapped genes. It ensures distinct mouse gene symbols are used.
h.genes <- m.genes %>%
  tibble::enframe("gene_index", "HGNC.symbol") %>%
  dplyr::left_join(human_to_mouse_homologs, by = "HGNC.symbol") %>%
  dplyr::distinct(HGNC.symbol, .keep_all = T)
h.genes <- h.genes[which(!(is.na(h.genes$MGI.symbol) | h.genes$MGI.symbol=="")),]
h.genes <- h.genes %>% distinct(MGI.symbol, .keep_all = T)

# This section creates the new Seurat object with mouse gene names
cd_human <- obj.h@assays$RNA@counts[h.genes$HGNC.symbol,]
rownames(cd_human) = h.genes$MGI.symbol
temp_meta <- obj.h@meta.data
obj.h$bc <- rownames(temp_meta)
temp_meta$nFeature_RNA <- NULL
temp_meta$nCount_RNA <- NULL
fobj.h <- CreateSeuratObject(cd_human)
fobj.h <- NormalizeData(fobj.h)
fobj.h@meta.data$bc <- rownames(fobj.h@meta.data)
temp_meta$bc <- rownames(temp_meta)
fobj.h@meta.data <- left_join(fobj.h@meta.data,temp_meta, by = "bc")
rownames(fobj.h@meta.data) <- fobj.h@meta.data$bc

saveRDS(fobj.h, "HumanMeninges_MouseGenes.rds")
