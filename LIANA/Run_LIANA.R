setwd("/n/scratch/users/p/pab1164/pab1164/ShamsLR/TGMeninges_v2")
library(Seurat)
library(liana)

# LIANA processing for Mouse

tg.mm <- readRDS("/n/scratch/users/p/pab1164/pab1164/ShamsLR/TG_Neurons_Mouse.rds")
men <- subset(skull, subset=region=="Meninges")
men <- subset(men, subset = level1 %in% unique(Idents(obj.h)))

Idents(men) <- men$level1
Idents(tg.mm) <- tg.hum$level2_final

DefaultAssay(men) <- "RNA"
DefaultAssay(tg.mm) <- "RNA"

obj <- merge(tg.mm, men, merge.data = TRUE)
table(Idents(obj))
obj <- subset(obj, downsample=1000)
gc()

obj.li.r <- liana_wrap(obj, method = c("logfc","sca"), resource = "MouseConsensus", return_all = TRUE, expr_prop=0.01, parallelize = TRUE, workers=8)
gc()
obj.agg.r <- liana_aggregate(obj.li.r, resource = "MouseConsensus")
saveRDS(obj.li.r, "TG_Meninges.Mouse.LianaWrapped_Relaxed.rds")
saveRDS(obj.agg.r, "TG_Meninges.Mouse.LianaAggregated_Relaxed.rds")

# LIANA processing for Human

obj.h <- readRDS("/n/scratch/users/p/pab1164/pab1164/ShamsLR/TGMeninges_v2/HumanMeninges.Anchored.Filtered.rds")
tg.hum <- readRDS("/n/scratch/users/p/pab1164/pab1164/ShamsLR/TG_Neurons_Humans.rds")

Idents(obj.h) <- obj.h$predicted.celltype
Idents(tg.hum) <- tg.hum$level2_final

DefaultAssay(obj.h) <- "RNA"
DefaultAssay(tg.hum) <- "RNA"

obj <- merge(tg.hum, obj.h, merge.data = TRUE)
obj <- subset(obj, downsample=1000)
gc()

obj.li.r <- liana_wrap(obj, method = c("logfc","sca"), resource = "MouseConsensus", return_all = TRUE, expr_prop=0.01, parallelize = TRUE, workers=8)
gc()
obj.agg.r <- liana_aggregate(obj.li.r, resource = "MouseConsensus")
saveRDS(obj.li.r, "TG_Meninges.Human.LianaWrapped_Relaxed.rds")
saveRDS(obj.agg.r, "TG_Meninges.Human.LianaAggregated_Relaxed.rds")
