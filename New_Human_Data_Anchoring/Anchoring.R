# This template script was used to anchor new human data to the atlas as described in the Methods section of the associated paper 

setwd("/mnt/data/projects/parth/Atlas_cellbend/MapCellBend/Seurat_pipelines/WashU/Shams_parameters/Anchored/Neurons/")
library(Seurat)

atlas <- readRDS("/mnt/data/projects/parth/Atlas_cellbend/MapCellBend/Seurat_pipelines/DRG.Neurons.Atlas.Primates.rds")
obj <- readRDS("/mnt/data/projects/parth/Atlas_cellbend/MapCellBend/Seurat_pipelines/WashU/Shams_parameters/Anchored/Neurons/WashU.Neurons.biomart.rds")

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures=2000)
obj <- ScaleData(obj, features=rownames(obj))
obj <- RunPCA(obj)

VariableFeatures(atlas@assays$RNA) <- VariableFeatures(atlas@assays$integrated)   # Here, we use the variable features from the integrated assay but use the RNA assay with raw counts for the purpose of RNA-to-RNA anchoring
DefaultAssay(atlas) = "RNA"

atlas <- ScaleData(atlas)
atlas <- RunPCA(atlas, reduction.name="pca.rna")

transfer.anchors <- FindTransferAnchors(reference = atlas, query = obj, reference.assay = "RNA", query.assay = "RNA", reduction="cca")
predictions <- TransferData(anchorset = transfer.anchors, refdata = atlas$final, dims=1:30, weight.reduction="cca")
obj <- AddMetaData(obj, metadata = predictions)
DefaultAssay(atlas) = "integrated"
atlas <- RunUMAP(atlas, dims = 1:30, reduction = "pca", return.model = TRUE)
obj <-  MapQuery(anchorset = transfer.anchors, reference = atlas, query = obj, refdata = list(celltype = "final"), reference.reduction = "pca", reduction.model = "umap", transferdata.args = list(weight.reduction="cca"), integrateembeddings.args = list(reductions="cca"))

saveRDS(obj, "WashU.Neurons.Coembedded.Anchored.FixedFinal.Primates.rds")
