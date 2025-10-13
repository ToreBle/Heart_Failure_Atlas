library(Seurat)
library(SeuratDisk)

Convert("./multinichenet/Heart_atlas_reduced.h5ad", "./multinichenet/Heart_atlas_reduced.h5Seurat", assay = "RNA", overwrite = T)
Single_data <- LoadH5Seurat("./multinichenet/Heart_atlas_reduced.h5Seurat", meta.data=FALSE, misc = FALSE)
meta_data = read.csv("./multinichenet/Heart_atlas_reduced_meta_data.csv", row.names=1)
Single_data <- AddMetaData(Single_data, metadata = meta_data)
Idents(Single_data) <- Single_data$Annotation_1
saveRDS(Single_data, "./multinichenet/Heart_atlas_reduced.rds")