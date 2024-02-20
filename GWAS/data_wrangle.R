library(MAGMA.Celltyping)
library(EWCE)
library(One2One)
library(Seurat)
library(tidyverse)
library(rlist)
library(R.utils)
library(foreach)
library(doParallel)
library(tm)


# load and format annotation
setwd("~")
seurat_mouse_neuron <- readRDS("seurat_relabelled.Rds")
seurat_mouse_nonneuron <- readRDS("Seurat.Rds")

## subset for mouse in TG
seurat_mouse_neuron <- subset(seurat_mouse_neuron, Species == "Mouse")
seurat_mouse_nonneuron <- subset(seurat_mouse_nonneuron, Species == "Mouse")

## have annotation in final column
seurat_mouse_neuron@meta.data$final <- seurat_mouse_neuron@meta.data$level2_final
seurat_mouse_nonneuron@meta.data$final <- seurat_mouse_nonneuron@meta.data$level2_round2

seurat_mouse <- merge(seurat_mouse_neuron,seurat_mouse_nonneuron)

## counts normalization 
counts_mouse <- seurat_mouse[["RNA"]]@counts
counts_mouse <- NormalizeData(counts_mouse, normalization.method = "RC", scale.factor = 10000)


## store subtype information in a list
meta_mouse <- seurat_mouse@meta.data
meta_mouse$subtype <- meta_mouse$final
meta_mouse <- meta_mouse %>% select(subtype)

meta_mouse$subtype <- gsub("-","_",meta_mouse$subtype)
meta_mouse$subtype <- gsub("\\+","_",meta_mouse$subtype)


annotLevel_mouse_subtype <- list(l1 = meta_mouse$subtype)

counts_mouse <- as.matrix(counts_mouse)

## drop non-1:1 orthologs, drops genes from an SCT expression matrix if they do not significantly vary between any celltypes. 
# Makes this decision based on use of an ANOVA (implemented with Limma).
count_mouse_Dropped_subtype <- drop.uninformative.genes(exp = counts_mouse, level2annot = meta_mouse$subtype)

ctd_mouse_subtype_path <- generate.celltype.data(count_mouse_Dropped_subtype,annotLevel_mouse_subtype,groupName = "mouse_DRG_neuron")

ctd_mouse_subtype_path <- load(file = ctd_mouse_subtype_path)
ctd_mouse_subtype <- ctd

## calculate the quantile groups for each celltype within the single cell dataset
ctd_mouse_subtype = prepare.quantile.groups(ctd_mouse_subtype,specificity_species="mouse",numberOfBins= 40)


dir.create(paste0(format(Sys.Date(), "%Y%m%d"),"_TG_neuron_nonneuron_mouse"))
setwd(paste0(format(Sys.Date(), "%Y%m%d"),"_TG_neuron_nonneuron_mouse"))

saveRDS(ctd_mouse_subtype,"gene_quantile_40.Rds")