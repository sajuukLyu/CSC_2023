
# MESSAGE -----------------------------------------------------------------
#
# author: Yulin Lyu, Cheng Li Lab, Peking University
# email: lvyulin@pku.edu.cn
# date: 2023.2.20
#
# ---

# load package ------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(data.table)

library(Seurat)

`%ni%` <- function(a, b) return(! a %in% b)

# WOT ---------------------------------------------------------------------

allSo <- readRDS("middata/hADSC_0618_scRNA/allSo_qc.rds")

# prepare input

usedSample <- c(
  "ADSC", "S1D0.5", "S1D2", "S1D4", "S1D8",
  "S2D4", "S2D8",
  "S3D0.33", "S3D0.67", "S3D1", "S3D2", "S3D4", "S3D6", "S3D8", "S3D12")

sampleTime <- c(
  0, 0.5, 2, 4, 8,
  8 + c(4, 8),
  8 + 8 + c(1/3, 2/3, 1, 2, 4, 6, 8, 12)
)

names(sampleTime) <- usedSample

allSo <- allSo[, allSo$sample %in% usedSample]

allSo$time <- sampleTime[allSo$sample]
table(allSo@meta.data[, c("time", "sample")])

wotExp <- allSo[["RNA"]]@data[VariableFeatures(allSo), ] %>% as.matrix %>% t
wotExp <- cbind(rownames(wotExp), wotExp)
colnames(wotExp)[1] <- "id"

write.table(wotExp, "middata/hADSC_scRNA/wot/wotExp.txt", sep = "\t", row.names = F, quote = F)

wotDay <- as.matrix(allSo$time)
rownames(wotDay) <- colnames(allSo)
wotDay <- cbind(rownames(wotDay), wotDay)
colnames(wotDay) <- c("id", "day")
write.table(wotDay, "middata/hADSC_scRNA/wot/wotDay.txt", sep = "\t", quote = F, row.names = F)

# perform WOT analysis according to developer instructions (https://broadinstitute.github.io/wot/)

# identify S3D12 celltype

end <- allSo[, allSo$sample == "S3D12"]
end %<>% FindVariableFeatures
end %<>% ScaleData %>% RunPCA
end %>% ElbowPlot(ndims = 50)
end %<>% FindNeighbors(dims = 1:20)
end %<>% FindClusters(resolution = 1)
end %<>% RunUMAP(dims = 1:20)

DimPlot(end, group.by = "seurat_clusters", label = T)

FeaturePlot(end, c("POU5F1", "SOX2", "NANOG", "KLF4"))
FeaturePlot(end, c("COL1A2", "DCN"))

end$celltype <- "other"
end$celltype[end$seurat_clusters %in% c()] <- "hCiPS"

# the celltypes were used to calculate right probability

allSo$CiPS <- readRDS("middata/hADSC_scRNA/CiPS_score.rds")[colnames(allSo)]

saveRDS(allSo, "middata/hADSC_scRNA/allSo_traj.rds")

# identify right trajectory -----------------------------------------------

allSo <- readRDS("middata/hADSC_scRNA/allSo_traj.rds")

allSo$traj <- "F"
allSo$traj[allSo$sample == "prime"] <- "prime"
allSo$traj[allSo$sample == "H1"] <- "ES"

# define right trajectory according to hCiPS score
allSo$traj[allSo$seurat_clusters %in% c()] <- "S0"
allSo$traj[allSo$seurat_clusters %in% c()] <- "S1"
allSo$traj[allSo$seurat_clusters %in% c()] <- "S2"
allSo$traj[allSo$seurat_clusters %in% c()] <- "S3"
allSo$traj[allSo$seurat_clusters %in% c()] <- "S4"
allSo$traj[allSo$seurat_clusters %in% c()] <- "S5"
allSo$traj[allSo$seurat_clusters %in% c()] <- "S6"
allSo$traj[allSo$seurat_clusters %in% c()] <- "S7"
allSo$traj[allSo$seurat_clusters %in% c()] <- "S8"
allSo$traj[allSo$seurat_clusters %in% c()] <- "S9"

saveRDS(allSo, "middata/hADSC_scRNA/allSo_traj.rds")

# integrate H1 ------------------------------------------------------------

# perform LSI projection to get LSI coordinates lsiData
# refer to https://github.com/GreenleafLab/MPAL-Single-Cell-2019

allSo[["lsi"]] <- CreateDimReducObject(
  lsiData[colnames(allSo), ], assay = "RNA", key = "LSI_"
)

umapRes <- uwot::umap(
  lsiData[, 1:20],
  n_neighbors = 40,
  min_dist = 0.3,
  metric = "euclidean",
  ret_model = T,
  verbose = T)

umapData <- umapRes$embedding
rownames(umapData) <- colnames(allSo)
colnames(umapData) <- str_c("UMAP_", 1:2)

allSo[["umap"]] <- CreateDimReducObject(
  umapData[colnames(allSo), ], assay = "RNA", key = "UMAP_"
)

saveRDS(allSo, "middata/hADSC_scRNA/allSo_traj.rds")


