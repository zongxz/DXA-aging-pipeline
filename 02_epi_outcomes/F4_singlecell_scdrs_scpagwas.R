# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 02_epi_outcomes
# 单细胞：scDRS 与 scPagwas 细胞类型关联分析
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(data.table)
# GSE169396 ====
library(Seurat)
library(Matrix)
library(dplyr)

base_dir <- "/path/to/epi/data/GSE169396"
samples <- c("GSM5201883_S1", "GSM5201884_S2", "GSM5201885_S3", "GSM5201886_S4")

seu_list <- list()
for (s in samples) {
  cat("Loading:", s, "\n")
  
  mat <- ReadMtx(
    mtx = file.path(base_dir, paste0(s, "_matrix.mtx.gz")),
    cells = file.path(base_dir, paste0(s, "_barcodes.tsv.gz")),
    features = file.path(base_dir, paste0(s, "_features.tsv.gz")),
    feature.column = 2
  )
  
  seu <- CreateSeuratObject(counts = mat, project = s, min.cells = 0, min.features = 0)
  seu$sample <- s
  seu_list[[s]] <- seu
  cat("genes =", nrow(seu), "cells =", ncol(seu), "\n")
}

sc <- merge(seu_list[[1]], y = seu_list[2:4], add.cell.ids = samples)
sc

sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")

VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "sample", ncol = 3, pt.size = 0)


library(dplyr)

# 2% / 98% quantiles
feature_q <- quantile(sc$nFeature_RNA,
                      probs = c(0.02, 0.98),
                      na.rm = TRUE)

count_q <- quantile(sc$nCount_RNA, probs = c(0.02, 0.98), na.rm = TRUE)

cat("nFeature 2%-98%:\n")
print(feature_q)

cat("\nnCount 2%-98%:\n")
print(count_q)


# Apply original QC
pass_qc <- sc$nFeature_RNA >= feature_q[1] & sc$nFeature_RNA <= feature_q[2] & 
  sc$nCount_RNA >= count_q[1] & sc$nCount_RNA <= count_q[2] & sc$percent.mt <= 10

cat("\nTotal cells:", ncol(sc), "\n")
cat("Passed QC :", sum(pass_qc), "\n")
cat("Removed   :", sum(!pass_qc), "\n")


# Per-sample retention
qc_result <- sc@meta.data %>%
  mutate(pass_qc = pass_qc) %>%
  group_by(sample) %>%
  summarise(
    total = n(),
    retained = sum(pass_qc),
    removed = sum(!pass_qc),
    retention = round(100 * mean(pass_qc), 2),
    .groups = "drop"
  )

print(qc_result)



# 原文 PCA
DefaultAssay(sc_qc) <- "RNA"
sc_qc <- NormalizeData(sc_qc,
                       normalization.method = "LogNormalize",
                       scale.factor = 10000)

sc_qc <- FindVariableFeatures(sc_qc, selection.method = "vst", nfeatures = 2000)
sc_qc <- ScaleData(sc_qc, features = VariableFeatures(sc_qc))
sc_qc <- RunPCA(sc_qc, features = VariableFeatures(sc_qc), npcs = 20)


sc_qc <- FindNeighbors(sc_qc, dims = 1:20)
sc_qc <- FindClusters(sc_qc, algorithm = 1)
sc_qc <- RunUMAP(sc_qc, dims = 1:2)


# PCA

# 过滤
lower_q <- quantile(sc$nCount_RNA, 0.02)
upper_q <- quantile(sc$nCount_RNA, 0.98)
sc_qc <- subset(sc,  subset = nCount_RNA > lower_q &  nCount_RNA < upper_q &  percent.mt < 10)
DefaultAssay(sc_qc) <- "RNA"

sc_qc <- NormalizeData(sc_qc,
                       normalization.method = "LogNormalize",
                       scale.factor = 10000)

sc_qc <- FindVariableFeatures(sc_qc, selection.method = "vst", nfeatures = 2000)
sc_qc <- ScaleData(sc_qc, features = VariableFeatures(sc_qc))
sc_qc <- RunPCA(sc_qc, features = VariableFeatures(sc_qc), npcs = 50)

ElbowPlot(sc_qc, ndims = 50)

# 查看累计方差
cumsum(Stdev(sc_qc, reduction = "pca")^2 / sum(Stdev(sc_qc, reduction = "pca")^2))

sc_qc <- FindNeighbors(sc_qc, dims = 1:20)
sc_qc <- FindClusters(sc_qc, resolution = 0.3)
sc_qc <- RunUMAP(sc_qc, dims = 1:20)

DimPlot(sc_qc, reduction = "umap", label = TRUE, repel = TRUE)

# Granulocyte 0 1 2
FeaturePlot(sc_qc, features = c('ITGAM', 'CEACAM8'), ncol = 2)
# T/NK cell 3
FeaturePlot(sc_qc, features = c('CD3D', 'CD3E', 'CD3G', 'NKG7', 'GNLY'), ncol = 3)
# mononuclear phagocyte  4 11
FeaturePlot(sc_qc, features = c('CD14', 'CD68'), ncol = 2)
# B cell 6 8 10 15
FeaturePlot(sc_qc, features = c('CD19', 'CD79A', 'MS4A1'), ncol = 3)
# Erythrocyte 7 5
FeaturePlot(sc_qc, features = c('GYPA', 'ALAS2'), ncol = 2)
# mesenchymal lin-eage cell 9
FeaturePlot( sc_qc, features = c('ALPL', 'DCN', 'COL1A1', 'THY1', 'SP7', 'DLX5', 'SOX9'), ncol = 4)
# Endothelial 13
FeaturePlot(sc_qc, features = c('VWF', 'PECAM1'), ncol = 2)
# Progenitor cell 12
FeaturePlot(sc_qc, features = c('KIT', 'PROM1'), ncol = 2)
# plasmacytoid dendritic cell 14
FeaturePlot(sc_qc, features = c('IL3RA', 'GZMB'), ncol = 2)


markers_check <- c('ITGAM', 'CEACAM8','CD3D', 'CD3E', 'CD3G', 'NKG7', 'GNLY','CD14', 
                   'CD68', 'CD19', 'CD79A', 'MS4A1','GYPA', 'ALAS2','ALPL', 'DCN', 
                   'COL1A1', 'THY1', 'SP7', 'DLX5', 'SOX9','VWF', 'PECAM1','KIT', 
                   'PROM1','IL3RA', 'GZMB')

DotPlot(sc_qc, features = markers_check, group.by = "seurat_clusters",dot.scale = 6) + RotatedAxis()

cluster_to_celltype <- c(
  "0" = "Granulocyte",
  "1" = "Granulocyte",
  "2" = "Granulocyte",
  "3" = "T_NK",
  "4" = "MP",
  "5" = "Erythrocyte",
  "6" = "B_cell",
  "7" = "Erythrocyte",
  "8" = "B_cell",
  "9" = "MesLin",  
  "10" = "B_cell",
  "11" = "MP",
  "12" = "Progenitor",
  "13" = "Endothelial",
  "14" = "pDC",
  "15" = "B_cell"
)

sc_qc$cell_type <- unname(cluster_to_celltype[as.character(sc_qc$seurat_clusters)])
Idents(sc_qc) <- sc_qc$cell_type
table(sc_qc$cell_type)

sc_qc$cell_type <- Idents(sc_qc)
table(sc_qc$cell_type)


DimPlot(sc_qc, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE) +
  ggtitle("Cell Type Annotation (Based on Article Markers)")

saveRDS(sc_qc, 'data/GSE169396_cell_type.rds')

sc_qc <- readRDS('data/GSE169396_cell_type.rds')

DimPlot(sc_qc, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE) 

library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)
sc_qc <- readRDS('data/GSE169396_cell_type.rds')


obj <- sc_qc

DefaultAssay(obj) <- "RNA"

# Seurat v5 多 layer 合并
obj <- JoinLayers(obj, assay = "RNA")
Layers(obj[["RNA"]])

# 转成 SingleCellExperiment
sce <- as.SingleCellExperiment(obj, assay = "RNA")

assayNames(sce)
dim(sce)

assay(sce, "counts")

cnt <- assay(sce, "counts")

summary(cnt@x)
all(cnt@x == round(cnt@x))

writeH5AD(sce, file = "data/GSE169396/GSE169396.h5ad", X_name = "counts")



# GSE284089空转=====
library(Seurat)
library(data.table)
library(Matrix)
library(dplyr)
library(ggplot2)
data_dir <- "/path/to/epi/data/GSE284089"

spatial <- Load10X_Spatial(
  data.dir =
    "/path/to/epi/data/GSE284089/GSM8677818_visium",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "GSM8677818"
)
head(spatial@meta.data)
spatial[["percent.mt"]] <- PercentageFeatureSet(spatial, pattern = "^MT-")

DefaultAssay(spatial) <- "Spatial"

cat("Number of spots:", ncol(spatial), "\n")
cat("Number of genes:", nrow(spatial), "\n")

cat("\nMetadata:\n")
head(spatial@meta.data)

cat("\nAssays:\n")
print(Assays(spatial))

cat("\nSpatial layers:\n")
print(Layers(spatial[["Spatial"]]))

summary(spatial$nCount_Spatial)

quantile(
  spatial$nCount_Spatial,
  probs = c(
    0, 0.01, 0.05, 0.10,
    0.25, 0.50, 0.75,
    0.90, 0.95, 0.99, 1
  )
)

summary(spatial$nFeature_Spatial)

quantile(
  spatial$nFeature_Spatial,
  probs = c(
    0, 0.01, 0.05, 0.10,
    0.25, 0.50, 0.75,
    0.90, 0.95, 0.99, 1
  )
)

summary(spatial$percent.mt)

quantile(
  spatial$percent.mt,
  probs = c(
    0, 0.25, 0.50, 0.75,
    0.90, 0.95, 0.99, 1
  )
)

spatial_qc <- subset(spatial, subset = nCount_Spatial >= 100)
cat("Before QC:", ncol(spatial), "\n")
cat("After QC :", ncol(spatial_qc), "\n")
cat("Removed   :", ncol(spatial) - ncol(spatial_qc), "\n")

saveRDS(spatial_qc, file = file.path(data_dir, "GSE284089_spatial_QC_raw.rds"))

# 转换h5ad
spatial <- readRDS("data/GSE284089/GSM8677818_spatial_QC_raw.rds")

DefaultAssay(spatial) <- "Spatial"

cat("Original dimension:\n")
print(dim(spatial))

# Check UMI QC
cat("\nUMI summary:\n")
print(summary(spatial$nCount_Spatial))

n_low <- sum(spatial$nCount_Spatial < 100)

cat("\nSpots with UMI < 100:", n_low, "\n")

cat("\nAfter QC:\n")
print(dim(spatial))

# Raw counts
counts <- LayerData(spatial, assay = "Spatial", layer = "counts")

cat("\nCounts dimension:\n")
print(dim(counts))

# Metadata
meta <- spatial@meta.data

# Spatial coordinates
coords <- GetTissueCoordinates(spatial)

cat("\nCoordinates:\n")
print(head(coords))
print(colnames(coords))

# Check barcode matching
stopifnot(all(colnames(spatial) == coords$cell))

# Create SCE
sce <- SingleCellExperiment(assays = list(counts = counts), colData = meta)

# Add spatial coordinates
reducedDim(sce, "spatial") <- as.matrix(coords[, c("x", "y")])

# Write h5ad
writeH5AD(sce, file = 'data/GSE284089/GSE284089_spatial.h5ad', X_name = "counts")




spatial <- readRDS("data/GSE284089/GSM8677818_spatial_QC_raw.rds")
# 对空间转录组数据进行 SCTransform 标准化
# assay = "Spatial" 表示使用空间转录组原始表达矩阵
# 标准化后的结果会保存到 SCT assay 中
spatial <- SCTransform(spatial, assay = "Spatial", verbose = FALSE)
# 基于 SCT assay 进行 PCA 降维# PCA 用于提取主要变异信息，为后续聚类和 UMAP 做准备
spatial <- RunPCA(spatial, assay = "SCT", verbose = FALSE)
# 根据 PCA 结果构建细胞/spot 之间的邻近关系
# dims = 1:30 表示使用前 30 个主成分
spatial <- FindNeighbors(spatial, reduction = "pca", dims = 1:30)
# 基于邻近图进行无监督聚类
# 聚类结果会保存到 Seurat 对象的 metadata 中
spatial <- FindClusters(spatial, verbose = FALSE)
# 基于 PCA 结果进行 UMAP 降维可视化# UMAP 可以把高维表达数据压缩到二维空间中展示
spatial <- RunUMAP(spatial, reduction = "pca", dims = 1:30)
spatial
# 在空间转录组切片上展示 ScType 的细胞类型注释结果
# group.by = "sctype_classification" 表示按照 ScType 注释结果着色
SpatialDimPlot(spatial, label = TRUE)

library(BayesSpace)
library(SingleCellExperiment)

# 从 Seurat 对象转换（需要 Spatial assay 和坐标）
sce <- as.SingleCellExperiment(spatial, assay = "Spatial")

# 预处理（Visium 平台）
sce <- spatialPreprocess(sce, platform = "Visium",
                         n.PCs = 15, n.HVGs = 2000, log.normalize = TRUE)

# 先用 qTune 选择 q 值
sce <- qTune(sce, qs = seq(3, 10), burn.in = 100, nrep = 1000)
qPlot(sce)   # 查看负对数似然曲线，选择拐点处的 q

# 运行空间聚类（示例 q=8, gamma=3 增强平滑）
sce <- spatialCluster(sce, q = 8, platform = "Visium",
                      d = 15, gamma = 3, nrep = 10000)

# 将结果放回 Seurat 对象以便可视化
spatial$BayesSpace_cluster <- sce$spatial.cluster
SpatialDimPlot(spatial, group.by = "BayesSpace_cluster", label = TRUE)








# GSE317068 ====
data <- readRDS(
  'data/GSE317068/GSE317068_Human_bone_scRNAseq_seurat_object_20251120.Rds'
)
DimPlot(data, reduction = "umap", group.by = 'Cluster_annotation')

data
colnames(data@meta.data)
unique(data@meta.data$Cluster_annotation)
head(rownames(data), 20)
head(Genes_by_pathway_kegg)
head(rownames(data))
# 检查基因名称是否与KEGG匹配
pathway_genes <- unique(unlist(Genes_by_pathway_kegg))

matched_genes <- intersect(pathway_genes, rownames(data))

length(pathway_genes)
length(matched_genes)
head(matched_genes)

Idents(data)
Idents(data) <- "Cluster_annotation"
Idents(data)


set.seed(20260903)

# 读取原始数据
data <- readRDS("data/GSE317068/GSE317068_Human_bone_scRNAseq_seurat_object_20251120.Rds")

# 设置细胞类型
Idents(data) <- "Cluster_annotation"

# 去除 Multiplets
data_clean <- subset(data, subset = Cluster_annotation != "18. Multiplets")

Idents(data_clean) <- "Cluster_annotation"
Idents(data_clean)
unique(data_clean@meta.data$Cluster_annotation)
saveRDS(data_clean, "data/GSE317068/GSE317068.rds")


data_clean$Celltype_scPagwas <- as.character(data_clean$Cluster_annotation)

# MSCs
data_clean$Celltype_scPagwas[data_clean$Cluster_annotation %in% c("0. MSCs CXCL12+GALNT17+", "1. MSCs CXCL12+STMN2+")] <- "MSCs"

# Osteoblasts
data_clean$Celltype_scPagwas[data_clean$Cluster_annotation %in% c("3. Pre-osteoblasts", "12. Mature osteoblasts")] <- "Osteoblasts"

# Chondrocytes / fibroblasts
data_clean$Celltype_scPagwas[data_clean$Cluster_annotation ==
                               "15. Chondrocytes and fibroblasts"] <- "Chondrocytes_Fibroblasts"

# Endothelial
data_clean$Celltype_scPagwas[data_clean$Cluster_annotation == "2. ECs"] <- "Endothelial"

# Pericytes
data_clean$Celltype_scPagwas[data_clean$Cluster_annotation == "6. Pericytes"] <- "Pericytes"

# SMCs
data_clean$Celltype_scPagwas[data_clean$Cluster_annotation == "7. SMCs"] <- "SMCs"

# Neutrophils
data_clean$Celltype_scPagwas[data_clean$Cluster_annotation %in% c(
  "4. Immature neutrophils",
  "8. Mature neutrophils IL1R2+",
  "17. Mature neutrophils LST1+"
)] <- "Neutrophils"

# Monocytes / macrophages
data_clean$Celltype_scPagwas[data_clean$Cluster_annotation ==
                               "11. Monocytes and macrophages"] <- "Monocytes_Macrophages"

# GMPs
data_clean$Celltype_scPagwas[data_clean$Cluster_annotation %in% c("10. GMPs DEFA1+", "14. GMPs CDC25C+")] <- "GMPs"

# B cells
data_clean$Celltype_scPagwas[data_clean$Cluster_annotation %in% c("9. Immature B cells", "16. Pre B cells", "19. B cells")] <- "B_cells"

# T / NK
data_clean$Celltype_scPagwas[data_clean$Cluster_annotation == "5. T and NK cells"] <- "T_NK"

# HSC
data_clean$Celltype_scPagwas[data_clean$Cluster_annotation == "13. HSCs"] <- "HSCs"

table(data_clean$Celltype_scPagwas)

DimPlot(data_clean, reduction = "umap", group.by = 'Celltype_scPagwas')

Idents(data_clean) <- "Celltype_scPagwas"


saveRDS(data_clean, "data/GSE317068/GSE317068_merge_cell.rds")


# 每个细胞类型随机抽取 10%
meta <- data_clean@meta.data

sample_cells <- meta %>%
  mutate(cell_id = rownames(meta)) %>%
  group_by(Cluster_annotation) %>%
  slice_sample(prop = 0.10) %>%
  pull(cell_id)

length(sample_cells)

# 构建测试数据
data_test <- data_clean[, sample_cells]

# 检查每个细胞类型数量
table(data_test$Cluster_annotation)


set.seed(20260903)

sample_cells <- meta %>%
  mutate(cell_id = rownames(meta)) %>%
  group_by(Cluster_annotation) %>%
  group_modify( ~ {
    n <- nrow(.x)
    n_sample <- max(30, ceiling(n * 0.10))
    n_sample <- min(n_sample, n)
    
    .x %>% slice_sample(n = n_sample)
  }) %>%
  pull(cell_id)

data_test <- data_clean[, sample_cells]

table(data_test$Cluster_annotation)
ncol(data_test)
sum(table(data_test$Cluster_annotation))
length(unique(data_test$Cluster_annotation))
saveRDS(data_test,
        "data/GSE317068/GSE317068_test_10percent.rds")



# GSE317068 scPagwas ====
load('result/F4_scPagwas/double_right/double_right_scPagwas.RData')
Pagwas$Celltype_scPagwas <- as.character(Pagwas$Cluster_annotation)
# MSCs
Pagwas$Celltype_scPagwas[Pagwas$Cluster_annotation %in% c("0. MSCs CXCL12+GALNT17+", "1. MSCs CXCL12+STMN2+")] <- "Mesenchymal stromal cells"
# Osteoblasts
Pagwas$Celltype_scPagwas[Pagwas$Cluster_annotation %in% c("3. Pre-osteoblasts", "12. Mature osteoblasts")] <- "Osteoblasts"
# Chondrocytes / fibroblasts
Pagwas$Celltype_scPagwas[Pagwas$Cluster_annotation == "15. Chondrocytes and fibroblasts"] <- "Chondrocytes fibroblasts"
# Endothelial
Pagwas$Celltype_scPagwas[Pagwas$Cluster_annotation == "2. ECs"] <- "Endothelial"
# Pericytes
Pagwas$Celltype_scPagwas[Pagwas$Cluster_annotation == "6. Pericytes"] <- "Pericytes"
# SMCs
Pagwas$Celltype_scPagwas[Pagwas$Cluster_annotation == "7. SMCs"] <- "Smooth muscle cells"
# Neutrophils
Pagwas$Celltype_scPagwas[Pagwas$Cluster_annotation %in% 
                           c("4. Immature neutrophils", "8. Mature neutrophils IL1R2+", 
                             "17. Mature neutrophils LST1+")] <- "Neutrophils"
# Monocytes / macrophages
Pagwas$Celltype_scPagwas[Pagwas$Cluster_annotation ==
                           "11. Monocytes and macrophages"] <- "Monocytes and macrophages"
# GMPs
Pagwas$Celltype_scPagwas[Pagwas$Cluster_annotation %in% 
                           c("10. GMPs DEFA1+", "14. GMPs CDC25C+")] <- "Granulocyte-monocyte progenitor"
# B cells
Pagwas$Celltype_scPagwas[Pagwas$Cluster_annotation %in% c("9. Immature B cells", 
                                                          "16. Pre B cells", "19. B cells")] <- "B cells"
# T / NK
Pagwas$Celltype_scPagwas[Pagwas$Cluster_annotation == "5. T and NK cells"] <- "T/NK cells"

# HSC
Pagwas$Celltype_scPagwas[Pagwas$Cluster_annotation == "13. HSCs"] <- "Hematopoietic stem cells"


## UMAP ====
library(Seurat)


# 提取umap数据
umap = Pagwas@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>%
  cbind(cell_type = Pagwas@meta.data$Celltype_scPagwas) %>%# 注释后的label信息 ，改为cell_type
  cbind(cell_type_full = Pagwas@meta.data$Celltype_scPagwas_full)
table(umap$cell_type)
table(umap$cell_type_full)

write.csv(umap,'result/F4_scPagwas/umap_data.csv', row.names = T)
umap <- fread('result/F4_scPagwas/umap_data.csv')
# 首先，计算UMAP坐标的边界值，以便在多个图层中重复使用
umap_min_x <- min(umap$UMAP_1)
umap_min_y <- min(umap$UMAP_2)


# 修改名称
umap$cell_type[umap$cell_type == "Monocytes macrophages"] <- "Monocytes/macrophages"
umap$cell_type[umap$cell_type == "Chondrocytes fibroblasts"] <- "Chondrocytes/fibroblasts"

# 查看修改后的结果
table(umap$cell_type)

# 计算每个类的UMAP坐标均值作为中心点
centers <- umap %>%
  group_by(cell_type) %>%
  dplyr::summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2), .groups = 'drop')


table(umap$cell_type)

allcolour <- c("#F13F71","#B6E580","#F1766D","#7a52ed","#F3B1A0","#FFD92F",
               "#17BECF","#42B540","#E59CC4FF","#Fd7F0c","#a0a003","#1F77B4FF" )

# 设置因子顺序
umap$cell_type <- factor(umap$cell_type)
centers$cell_type <- factor(centers$cell_type)


# 初始化ggplot对象
p <- ggplot(umap, aes(x = UMAP_1, y = UMAP_2, color = cell_type_full)) +
  geom_point(size = 0.1, alpha = 1) +
  scale_color_manual(values = allcolour) +
  # 在每个类的中心添加标号
  geom_text(data = centers, aes(x = UMAP_1, y = UMAP_2, label = as.character(cell_type)),
            size = 6, color = 'black',
            hjust = 0.5, vjust = 0.5, family = "Times New Roman", fontface='bold') +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = 'white'),
    plot.background = element_rect(fill = "white"),
    legend.position = 'none'
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))  +
  geom_segment(
    aes(x = umap_min_x-2, y = umap_min_y, xend = umap_min_x + 3, yend = umap_min_y),
    colour = "black", linewidth = 0.6, arrow = arrow(length = unit(0.3, "cm"))
  ) +
  geom_segment(
    aes(x = umap_min_x-2, y = umap_min_y, xend = umap_min_x -2, yend = umap_min_y + 6),
    colour = "black", linewidth = 0.6, arrow = arrow(length = unit(0.3, "cm"))
  ) +
  annotate("text", x = umap_min_x + 0.5, y = umap_min_y - 1, label = "UMAP1",
           color = "black", size = 6,family = "Times New Roman", fontface='bold') +
  annotate("text", x = umap_min_x - 2.5, y = umap_min_y + 2.5, label = "UMAP2",
           color = "black", size = 6, angle = 90,family = "Times New Roman",fontface='bold')
p

ggsave('result/F4_scPagwas/cell_umap.png',p, width=9.5,height=9, dpi = 600)


# 计算每个类的UMAP坐标均值作为中心点
centers <- umap %>%
  group_by(cell_type_full) %>%
  dplyr::summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2), .groups = 'drop')


table(umap$cell_type_full)

allcolour <- c("#F13F71","#B6E580","#F1766D","#7a52ed","#F3B1A0","#FFD92F",
               "#17BECF","#42B540","#E59CC4FF","#Fd7F0c","#a0a003","#1F77B4FF" )

# 设置因子顺序
umap$cell_type_full <- factor(umap$cell_type_full)
centers$cell_type_full <- factor(centers$cell_type_full)


# 初始化ggplot对象
p <- ggplot(umap, aes(x = UMAP_1, y = UMAP_2, color = cell_type_full)) +
  geom_point(size = 0.1, alpha = 1) +
  scale_color_manual(values = allcolour) +
  # 在每个类的中心添加标号
  geom_text(data = centers, aes(x = UMAP_1, y = UMAP_2, label = as.character(cell_type_full)),
            size = 4, color = 'black',
            hjust = 0.5, vjust = 0.5, family = "Times New Roman", fontface='bold') +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = 'white'),
    plot.background = element_rect(fill = "white"),
    legend.position = 'none'
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))  +
  geom_segment(
    aes(x = umap_min_x-2, y = umap_min_y, xend = umap_min_x + 3, yend = umap_min_y),
    colour = "black", linewidth = 0.6, arrow = arrow(length = unit(0.3, "cm"))
  ) +
  geom_segment(
    aes(x = umap_min_x-2, y = umap_min_y, xend = umap_min_x -2, yend = umap_min_y + 6),
    colour = "black", linewidth = 0.6, arrow = arrow(length = unit(0.3, "cm"))
  ) +
  annotate("text", x = umap_min_x + 0.5, y = umap_min_y - 1, label = "UMAP1",
           color = "black", size = 4,family = "Times New Roman", fontface='bold') +
  annotate("text", x = umap_min_x - 2.5, y = umap_min_y + 2.5, label = "UMAP2",
           color = "black", size = 4, angle = 90,family = "Times New Roman",fontface='bold')
p

ggsave('result/F4_scPagwas/cell_umap_full.png',p, width=8,height=7, dpi = 600)

## scPagwas小提琴图====
load('result/F4_scPagwas/double_right/double_right_scPagwas.RData')
double_meta <- Pagwas@meta.data
double_meta$cell <- rownames(double_meta)
write.csv(double_meta, 'result/F4_scPagwas/double_right_meta.csv',row.names = F)

load('result/F4_scPagwas/single_right_hip/single_right_hip_scPagwas.RData')
hip_meta <- Pagwas@meta.data
hip_meta$cell <- rownames(hip_meta)
write.csv(hip_meta, 'result/F4_scPagwas/hip_meta.csv',row.names = F)

load('result/F4_scPagwas/single_right_knee/single_right_knee_scPagwas.RData')
knee_meta <- Pagwas@meta.data
knee_meta$cell <- rownames(knee_meta)
write.csv(knee_meta, 'result/F4_scPagwas/knee_meta.csv',row.names = F)



# meta data
double_meta <- fread('result/F4_scPagwas/double_right_meta.csv')
knee_meta <- fread('result/F4_scPagwas/knee_meta.csv')
hip_meta <- fread('result/F4_scPagwas/hip_meta.csv')
names(double_meta)
names(knee_meta)
names(hip_meta)
# P value
double_pval <- fread('result/F4_scPagwas/double_right/double_right_Merged_celltype_pvalue.csv')
knee_pval <- fread('result/F4_scPagwas/single_right_knee/single_right_knee_Merged_celltype_pvalue.csv')
hip_pval <- fread('result/F4_scPagwas/single_right_hip/single_right_hip_Merged_celltype_pvalue.csv')
names(double_pval)
names(knee_pval)
names(hip_pval)


double_pval$celltype
table(double_meta$Celltype_scPagwas)
table(double_meta$Celltype_scPagwas_full)

allcolour <- c("#F13F71","#B6E580","#F1766D","#7a52ed","#F3B1A0","#FFD92F",
               "#17BECF","#42B540","#E59CC4FF","#Fd7F0c","#a0a003","#1F77B4FF" )
map_df <- data.frame(
  celltype            = c("B_cells", "Chondrocytes_Fibroblasts", "Endothelial",
                          "GMPs", "HSCs", "Monocytes_Macrophages", "MSCs",
                          "Neutrophils", "Osteoblasts", "Pericytes", "SMCs", "T_NK"),
  Celltype_scPagwas      = c("B cells", "Chondrocytes/fibroblasts", "ECs",
                             "GMPs", "HSCs", "Monocytes/macrophages", "MSCs",
                             "Neutrophils", "Osteoblasts", "Pericytes", "SMCs", "T/NK cells"),
  Celltype_scPagwas_full = c("B cells", "Chondrocytes/fibroblasts", "Endothelial",
                             "Granulocyte-monocyte progenitor", "Hematopoietic stem cells",
                             "Monocytes and macrophages", "Mesenchymal stromal cells",
                             "Neutrophils", "Osteoblasts", "Pericytes",
                             "Smooth muscle cells", "T/NK cells")
)

double_pval <- double_pval %>%
  left_join(map_df, by = "celltype")
knee_pval <- knee_pval %>%
  left_join(map_df, by = "celltype")
hip_pval <- hip_pval %>%
  left_join(map_df, by = "celltype")


knee_pval
table(knee_meta$Celltype_scPagwas)
# meta data
double_meta <- fread('result/F4_scPagwas/double_right_meta.csv')
knee_meta <- fread('result/F4_scPagwas/knee_meta.csv')
hip_meta <- fread('result/F4_scPagwas/hip_meta.csv')
names(double_meta)
names(knee_meta)
names(hip_meta)
# P value
double_pval <- fread('result/F4_scPagwas/double_right/double_right_Merged_celltype_pvalue.csv')
knee_pval <- fread('result/F4_scPagwas/single_right_knee/single_right_knee_Merged_celltype_pvalue.csv')
hip_pval <- fread('result/F4_scPagwas/single_right_hip/single_right_hip_Merged_celltype_pvalue.csv')

double_pval$celltype
table(double_meta$Celltype_scPagwas)
table(double_meta$Celltype_scPagwas_full)
allcolour <- c("#F13F71","#B6E580","#F1766D","#7a52ed","#F3B1A0","#FFD92F",
               "#17BECF","#42B540","#E59CC4FF","#Fd7F0c","#a0a003","#1F77B4FF" )

map_df <- data.frame(
  celltype            = c("B_cells", "Chondrocytes_Fibroblasts", "Endothelial",
                          "GMPs", "HSCs", "Monocytes_Macrophages", "MSCs",
                          "Neutrophils", "Osteoblasts", "Pericytes", "SMCs", "T_NK"),
  Celltype_scPagwas      = c("B cells", "Chondrocytes/fibroblasts", "ECs",
                             "GMPs", "HSCs", "Monocytes/macrophages", "MSCs",
                             "Neutrophils", "Osteoblasts", "Pericytes", "SMCs", "T/NK cells"),
  Celltype_scPagwas_full = c("B cells", "Chondrocytes fibroblasts", "Endothelial",
                             "Granulocyte-monocyte progenitor", "Hematopoietic stem cells",
                             "Monocytes and macrophages", "Mesenchymal stromal cells",
                             "Neutrophils", "Osteoblasts", "Pericytes",
                             "Smooth muscle cells", "T/NK cells")
)

double_pval <- double_pval %>%
  left_join(map_df, by = "celltype")
knee_pval <- knee_pval %>%
  left_join(map_df, by = "celltype")
hip_pval <- hip_pval %>%
  left_join(map_df, by = "celltype")


knee_pval
table(knee_meta$Celltype_scPagwas)


# Knee
# 生成 p 值显示文本
# 统一设置 P 值高度
y_p <- max(knee_meta$scPagwas.TRS.Score1, na.rm = TRUE) +
  0.08 * diff(range(knee_meta$scPagwas.TRS.Score1, na.rm = TRUE))

knee_pval_label <- knee_pval %>%
  mutate(p_label = ifelse(
    pvalue < 0.001, "bolditalic(P)~bold(\" < 0.001\")",
    paste0("bolditalic(P)~bold(\" = ", sprintf('%.3f', pvalue), "\")")),
    y = y_p) %>%
  dplyr::select(Celltype_scPagwas, p_label, pvalue, y)


knee_meta$trait <- knee_pval_label$trait <- 'Knee'


# Hip
y_p <- max(hip_meta$scPagwas.TRS.Score1, na.rm = TRUE) +
  0.08 * diff(range(hip_meta$scPagwas.TRS.Score1, na.rm = TRUE))

hip_pval_label <- hip_pval %>%
  mutate(p_label = ifelse(
    pvalue < 0.001, "bolditalic(P)~bold(\" < 0.001\")",
    paste0("bolditalic(P)~bold(\" = ", sprintf('%.3f', pvalue), "\")")),
    y = y_p) %>%
  dplyr::select(Celltype_scPagwas, p_label, pvalue, y)

hip_meta$trait <- hip_pval_label$trait <- 'Hip'


# Knee-Hip
y_p <- max(double_meta$scPagwas.TRS.Score1, na.rm = TRUE) +
  0.08 * diff(range(double_meta$scPagwas.TRS.Score1, na.rm = TRUE))

double_pval_label <- double_pval %>%
  mutate(p_label = ifelse(
    pvalue < 0.001, "bolditalic(P)~bold(\" < 0.001\")",
    paste0("bolditalic(P)~bold(\" = ", sprintf('%.3f', pvalue), "\")")),
    y = y_p) %>%
  dplyr::select(Celltype_scPagwas, p_label, pvalue, y)

double_meta$trait <- double_pval_label$trait <- 'Knee - Hip'


meta <- rbind(knee_meta, hip_meta)
meta <- rbind(meta, double_meta)


pval_label <- rbind(knee_pval_label, hip_pval_label)
pval_label <- rbind(pval_label, double_pval_label)
pval_label$y <- max(pval_label$y) - 0.04

names(meta)
meta <- meta[,c("scPagwas.TRS.Score1","Celltype_scPagwas", "Celltype_scPagwas_full",  "trait")]

# 修改名称
meta$Celltype_scPagwas[meta$Celltype_scPagwas == "Monocytes macrophages"] <- "Monocytes/macrophages"
meta$Celltype_scPagwas[meta$Celltype_scPagwas == "Chondrocytes fibroblasts"] <- "Chondrocytes/fibroblasts"

head(meta)
head(pval_label)
meta$trait <- factor(meta$trait, levels = c('Knee', 'Hip', 'Knee - Hip'))
pval_label$trait <- factor(pval_label$trait, levels = c('Knee', 'Hip', 'Knee - Hip'))


trs_median <- meta %>%
  group_by(trait, Celltype_scPagwas,Celltype_scPagwas_full) %>%
  summarise(
    median_TRS = median(scPagwas.TRS.Score1, na.rm = TRUE),
    n_cells = sum(!is.na(scPagwas.TRS.Score1)),
    .groups = "drop"
  ) %>%
  arrange(trait, Celltype_scPagwas,Celltype_scPagwas_full)

print(trs_median)

meta_pval <- inner_join(trs_median, pval_label)
write.csv(meta_pval, 'result/F4_scPagwas/cell_trs_pval.csv')


p <- ggplot(meta, aes(x = Celltype_scPagwas, y = scPagwas.TRS.Score1, fill = Celltype_scPagwas)) +
  geom_violin(scale = "width", width = 0.5, trim = TRUE, linewidth = 0.6, color = "black") +
  geom_boxplot(width = 0.09, outlier.shape = NA, linewidth = 0.5, color = "black", fill = "white") +
  stat_summary(fun = median, geom = "point", size = 0.8, color = "black") +
  geom_text(data = pval_label,
            aes(x = Celltype_scPagwas, y = y, label = p_label),
            inherit.aes = FALSE, parse = TRUE,
            size = 4.5, family = "Times New Roman", colour = "black") +
  scale_fill_manual(values = allcolour) +
  scale_y_continuous(limits = c(min(meta$scPagwas.TRS.Score1), 1.224)) +
  facet_wrap(~ trait, nrow = 3, scales = "fixed") +
  labs(x = NULL, y = "Trait-relevant score", fill = NULL) +
  theme_classic() +
  theme(
    panel.grid   = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.background = element_rect(fill = "grey85", color = "black", linewidth = 1),
    strip.text   = element_text(size = 16, family = "Times New Roman", face = "bold", color = "black"),
    axis.title.y = element_text(size = 20, family = "Times New Roman", face = "bold", color = "black"),
    axis.text.y  = element_text(size = 18, family = "Times New Roman", face = "bold", color = "black"),
    axis.title.x = element_blank(),
    panel.spacing.y  = unit(0.05, "cm"),
    axis.text.x  = element_text(size = 16, family = "Times New Roman", face = "bold", color = "black", 
                                angle = 40, hjust = 1),
    legend.position = "none"
  )
p
ggsave('result/F4_scPagwas/volin_cell.png',p, width=12, height=9, dpi = 600)

# GSE317068 scDRS ====

## MAGMA 结果处理 ====
root_dir <- "result/F3_FUMA"
sub_dirs <- list.dirs(path = root_dir, full.names = TRUE, recursive = FALSE)

for (dir in sub_dirs) {
  gene_res_path <- file.path(dir, "POS", "magma.genes.out")
  
  if (file.exists(gene_res_path)) {
    
    gene_res <- read.table(gene_res_path, header = TRUE)
    zscore_file <- gene_res %>%
      transmute(
        GENE = SYMBOL,
        sleepscore = ZSTAT
      ) %>% 
      distinct(GENE, .keep_all = TRUE)
    
    # 构建输出文件路径
    output_path <- file.path("result/F4_scDRS/data",
                             paste0(basename(dir), "_trait_gene_zscore.tsv"))
    
    # 写入文件
    write.table(zscore_file, output_path, sep = "\t", quote = FALSE, row.names = FALSE)
    
    print(paste("Processed", dir))
  } else {
    print(paste("File not found in", dir))
  }
}

## 转h5文件 ====

library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)


data_test <- readRDS("data/GSE317068/GSE317068_merge_cell.rds")
# 查看对象
class(data_test)
data_test
Assays(data_test)
DefaultAssay(data_test)

data_test[["RNA"]]
Layers(data_test[["RNA"]])


# Seurat -> SingleCellExperiment
sce <- as.SingleCellExperiment(data_test, assay = "RNA")

# 查看转换后的 assay
assayNames(sce)
dim(sce)

writeH5AD(sce, file = "data/GSE317068/GSE317068_merge_cell.h5ad", X_name = "counts")


## scDRS 结果 ====
library(data.table)
# scDRS umap
score <- fread('/path/to/epi/result/F4_scDRS/result/double.score.gz')
fullscore <- fread('/path/to/epi/result/F4_scDRS/result/double.full_score.gz')
scdrs_group <- fread('/path/to/epi/result/F4_scDRS/result/double.scdrs_group.Celltype_scPagwas')
scdrs_gene <- fread('/path/to/epi/result/F4_scDRS/result/double.scdrs_gene')




double <- fread('/path/to/epi/result/F4_scDRS/result/double.scdrs_group.Celltype_scPagwas')
hip <- fread('/path/to/epi/result/F4_scDRS/result/hip.scdrs_group.Celltype_scPagwas')
knee <- fread('/path/to/epi/result/F4_scDRS/result/knee.scdrs_group.Celltype_scPagwas')





meta <- Pagwas@meta.data
meta$V1 <- row.names(meta)
colnames(meta)

scdrs_score <- inner_join(score, meta[,c(17, 10)])
colnames(scdrs_score)

# 取 UMAP 坐标
umap_df <- as.data.frame(Pagwas@reductions$umap@cell.embeddings)
umap_df$V1 <- rownames(umap_df)

colnames(umap_df) <- c("UMAP_1", "UMAP_2", "V1")

# 合并 scDRS 结果
plot_df <- left_join(umap_df, scdrs_score, by = "V1")

#  绘制 UMAP（颜色 = norm_score）
ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = norm_score)) +
  geom_point(size = 0.6, alpha = 0.8) +
  scale_color_gradient(low = '#fce0d4', high = 'red') +
  theme_minimal() +
  labs(color = "norm_score")


