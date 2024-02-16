# This is the template script used for any anchoring done between objects for the purpose of LR-interactions using LIANA only.  

setwd("~/scratch_pab/pab1164/Shams_jobs/Skin_DRG_LR")
library(Seurat)

obj.m <- readRDS("Skin.Mouse.Prepared.Object.rds")
obj.h <- readRDS("Skin.Human.Prepared.Object.rds")
obj.m <- FindVariableFeatures(obj.m,selection.method = "vst")
obj.h <- FindVariableFeatures(obj.h,selection.method = "vst")
obj.m <- ScaleData(obj.m, features = rownames(obj.m))
obj.h <- ScaleData(obj.h, features = rownames(obj.h))
obj.m <- RunPCA(obj.m)
obj.h <- RunPCA(obj.h)
obj.m <- FindNeighbors(obj.m, dims = 1:30)
obj.h <- FindNeighbors(obj.h, dims=1:30)

anchors <-  FindTransferAnchors(reference = obj.h, query = obj.m, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = obj.h$cell_type, dims = 1:30)
obj.m <- AddMetaData(obj.m, metadata = predictions)
obj.h <- RunUMAP(object = obj.h, dims = 1:30, reduction = "pca", return.model = TRUE)
obj.m <- MapQuery(anchorset = anchors, reference = obj.h, query = obj.m,
                            refdata = list(celltype = "cell_type"), reference.reduction = "pca", reduction.model = "umap")
obj.m$pred_pass <- obj.m$prediction.score.max >= 0.5
sub <- subset(obj.m, subset = pred_pass==TRUE)
saveRDS(obj.m, "Mouse.Skin.Anchored.FullObj.rds")
saveRDS(sub, "Mouse.Skin.Anchored.Filtered.rds")
