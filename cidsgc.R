library(tidyverse)
library(Seurat)
library(ggpubr)
library(ggvenn)
library(ComplexHeatmap)
library(ProjecTILs)

# Read data ----
raw.data <- Read10X(data.dir = '/db2/users/shnam/04.IMP/10.CIDSGC/matrix_files/')
metadata <- read.csv("/db2/users/shnam/04.IMP/10.CIDSGC/metadata.csv")

# Data preprocessing ----
# X is a column for cell_id
rownames(metadata) <- metadata$X
metadata <- metadata[,-1]

# QC for meta data ----
metadata %>% 
  count(state) %>%
  ggplot(aes(reorder(state, (-n)), n)) + geom_bar(stat = "identity") + xlab("state") + ylab("# of cell")
metadata %>% 
  count(condition) %>% 
  ggplot(aes(reorder(condition, (-n)), n)) + geom_bar(stat = "identity") + xlab("condition") + ylab("# of cell") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Create Seurat object ----
sobj <- CreateSeuratObject(counts = raw.data, meta.data = metadata)

# Calculate percentage of reads that map to the mitochondrial genome ----
sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^mt-")

# Split SeuratObject by condition
sobj.list <- SplitObject(sobj, split.by = "condition")
sobj.list$Unperturbed

# Visualize QC metrics as a violin plot ----
VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(sobj.list$Unperturbed, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(sobj.list$Tox2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(sobj.list$Arid5b, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(sobj.list$Dvl2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Visualize feature-feature relationships ----
p1 <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1 + p2
p1 <- FeatureScatter(sobj.list$Unperturbed, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(sobj.list$Unperturbed, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1 + p2
p1 <- FeatureScatter(sobj.list$Tox2, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(sobj.list$Tox2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1 + p2
p1 <- FeatureScatter(sobj.list$Arid5b, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(sobj.list$Arid5b, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1 + p2
p1 <- FeatureScatter(sobj.list$Dvl2, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(sobj.list$Dvl2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1 + p2

# Normalizing the data ----
sobj <- NormalizeData(sobj)
sobj.unperturbed <- NormalizeData(sobj.list$Unperturbed)
sobj.tox2 <- NormalizeData(sobj.list$Tox2)
sobj.arid5b <- NormalizeData(sobj.list$Arid5b)
sobj.dvl2 <- NormalizeData(sobj.list$Dvl2)

# Identification of highly variable features ----
sobj <- FindVariableFeatures(sobj, nfeatures = 2000)
top10 <- head(VariableFeatures(sobj), 10)
p0 <- LabelPoints(plot = VariableFeaturePlot(sobj), points = top10, repel = T, xnudge = 0, ynudge = 0)
sobj.unperturbed <- FindVariableFeatures(sobj.unperturbed, nfeatures = 2000)
top10.unperturbed <- head(VariableFeatures(sobj.unperturbed), 10)
p1 <- LabelPoints(plot = VariableFeaturePlot(sobj.unperturbed), points = top10.unperturbed, repel = T, xnudge = 0, ynudge = 0)
sobj.tox2 <- FindVariableFeatures(sobj.tox2, nfeatures = 2000)
top10.tox2 <- head(VariableFeatures(sobj.tox2), 10)
p2 <- LabelPoints(plot = VariableFeaturePlot(sobj.tox2), points = top10.tox2, repel = T, xnudge = 0, ynudge = 0)
sobj.arid5b <- FindVariableFeatures(sobj.arid5b, nfeatures = 2000)
top10.arid5b <- head(VariableFeatures(sobj.arid5b), 10)
p3 <- LabelPoints(plot = VariableFeaturePlot(sobj.arid5b), points = top10.arid5b, repel = T, xnudge = 0, ynudge = 0)
sobj.dvl2 <- FindVariableFeatures(sobj.dvl2, nfeatures = 2000)
top10.dvl2 <- head(VariableFeatures(sobj.dvl2), 10)
p4 <- LabelPoints(plot = VariableFeaturePlot(sobj.dvl2), points = top10.dvl2, repel = T, xnudge = 0, ynudge = 0)
(p1+p2+p3+p4) & NoLegend()

top10 <- head(VariableFeatures(sobj), 100)
top10.unperturbed <- head(VariableFeatures(sobj.unperturbed), 100)
top10.tox2 <- head(VariableFeatures(sobj.tox2), 100)
top10.arid5b <- head(VariableFeatures(sobj.arid5b), 100)
top10.dvl2 <- head(VariableFeatures(sobj.dvl2), 100)
top10.list <- list(`Unperturbed` = top10.unperturbed,
                   `Tox2` = top10.tox2,
                   `Arid5b` = top10.arid5b,
                   `Dvl2` = top10.dvl2)
ggvenn(top10.list, c("Unperturbed", "Tox2", "Arid5b", "Dvl2"), show_percentage = F)

# Scaling the data ----
genes <- rownames(sobj)
sobj <- ScaleData(sobj, features = genes)
genes.unperturbed <- rownames(sobj.unperturbed)
sobj.unperturbed <- ScaleData(sobj.unperturbed, features = genes.unperturbed)
genes.tox2 <- rownames(sobj.tox2)
sobj.tox2 <- ScaleData(sobj.tox2, features = genes.tox2)
genes.arid5b <- rownames(sobj.arid5b)
sobj.arid5b <- ScaleData(sobj.arid5b, features = genes.arid5b)
genes.dvl2 <- rownames(sobj.dvl2)
sobj.dvl2 <- ScaleData(sobj.dvl2, features = genes.dvl2)

# Perform linear dimensional reduction ----
sobj <- RunPCA(sobj, features = VariableFeatures(object = sobj))

sobj.unperturbed <- RunPCA(sobj.unperturbed, features = VariableFeatures(object = sobj.unperturbed))
VizDimLoadings(sobj.unperturbed, dims = 1:2, reduction = "pca")
DimPlot(sobj.unperturbed, reduction = "pca")
DimHeatmap(sobj.unperturbed, dims = 1, cells = 500, balanced = T)
DimHeatmap(sobj.unperturbed, dims = 1:30, cells = 500, balanced = T)

sobj.tox2 <- RunPCA(sobj.tox2, features = VariableFeatures(object = sobj.tox2))

sobj.arid5b <- RunPCA(sobj.arid5b, features = VariableFeatures(object = sobj.arid5b))

sobj.dvl2 <- RunPCA(sobj.dvl2, features = VariableFeatures(object = sobj.dvl2))

# Determine the dimensionality of the dataset ----
sobj <- JackStraw(sobj, num.replicate = 100)
sobj <- ScoreJackStraw(sobj, dims = 1:20)
JackStrawPlot(sobj, dims = 1:20)
ElbowPlot(sobj)

sobj.unperturbed <- JackStraw(sobj.unperturbed, num.replicate = 100)
sobj.unperturbed <- ScoreJackStraw(sobj.unperturbed, dims = 1:20)
JackStrawPlot(sobj.unperturbed, dims = 1:20)
ElbowPlot(sobj.unperturbed)

sobj.tox2 <- JackStraw(sobj.tox2, num.replicate = 100)
sobj.tox2 <- ScoreJackStraw(sobj.tox2, dims = 1:20)
JackStrawPlot(sobj.tox2, dims = 1:20)
ElbowPlot(sobj.tox2)

sobj.arid5b <- JackStraw(sobj.arid5b, num.replicate = 100)
sobj.arid5b <- ScoreJackStraw(sobj.arid5b, dims = 1:20)
JackStrawPlot(sobj.arid5b, dims = 1:20)
ElbowPlot(sobj.arid5b)

sobj.dvl2 <- JackStraw(sobj.dvl2, num.replicate = 100)
sobj.dvl2 <- ScoreJackStraw(sobj.dvl2, dims = 1:20)
JackStrawPlot(sobj.dvl2, dims = 1:20)
ElbowPlot(sobj.dvl2)

# Cluster the cells ----
sobj <- FindNeighbors(sobj, dims = 1:15)
sobj <- FindClusters(sobj, resolution = 0.5)

sobj.unperturbed <- FindNeighbors(sobj.unperturbed, dims = 1:15)
sobj.unperturbed <- FindClusters(sobj.unperturbed, resolution = 0.5)

sobj.tox2 <- FindNeighbors(sobj.tox2, dims = 1:15)
sobj.tox2 <- FindClusters(sobj.tox2, resolution = 0.5)

sobj.arid5b <- FindNeighbors(sobj.arid5b, dims = 1:15)
sobj.arid5b <- FindClusters(sobj.arid5b, resolution = 0.5)

sobj.dvl2 <- FindNeighbors(sobj.dvl2, dims = 1:15)
sobj.dvl2 <- FindClusters(sobj.dvl2, resolution = 0.5)

# Run non-linear dimentional reduction ----
sobj <- RunUMAP(sobj, dims = 1:15)
p1 <- DimPlot(sobj, reduction = "umap", label = T) + ggtitle("All")
p2 <- DimPlot(sobj, reduction = "umap", label = T, group.by = "state")
DimPlot(sobj, reduction = "umap", label = T, group.by = "condition") & NoLegend()
DimPlot(sobj, reduction = "umap", label = T, group.by = "lane")
p1|p2

sobj.unperturbed <- RunUMAP(sobj.unperturbed, dims = 1:15)
p1 <- DimPlot(sobj.unperturbed, reduction = "umap", label = T) + ggtitle("Unperturbed")
p2 <- DimPlot(sobj.unperturbed, reduction = "umap", label = T, group.by = "state")
p1|p2

sobj.tox2 <- RunUMAP(sobj.tox2, dims = 1:15)
p1 <- DimPlot(sobj.tox2, reduction = "umap", label = T) + ggtitle("Tox2")
p2 <- DimPlot(sobj.tox2, reduction = "umap", label = T, group.by = "state")
p1|p2

sobj.arid5b <- RunUMAP(sobj.arid5b, dims = 1:15)
p1 <- DimPlot(sobj.arid5b, reduction = "umap", label = T) + ggtitle("Arid5b")
p2 <- DimPlot(sobj.arid5b, reduction = "umap", label = T, group.by = "state")
p1|p2

sobj.dvl2 <- RunUMAP(sobj.dvl2, dims = 1:15)
p1 <- DimPlot(sobj.dvl2, reduction = "umap", label = T) + ggtitle("Dvl1")
p2 <- DimPlot(sobj.dvl2, reduction = "umap", label = T, group.by = "state")
p1|p2

# Finding differentially expressed features ----
# Find markers for every cluster compared to all remaining cells

# There are several marker genes that are commonly used to identify T-cell specific expression in mice.
# Cd3e, Cd4, Cd8a, Tcrb (X), Lck
p1 <- VlnPlot(sobj.unperturbed, features = c("Cd3e")) + ggtitle("Unperturbed - Cd3e")
p2 <- VlnPlot(sobj.tox2, features = c("Cd3e")) + ggtitle("Tox2 - Cd3e")
p3 <- VlnPlot(sobj.arid5b, features = c("Cd3e")) + ggtitle("Arid5b - Cd3e")
p4 <- VlnPlot(sobj.dvl2, features = c("Cd3e")) + ggtitle("Dvl2 - Cd3e")
(p1|p2|p3|p4) & NoLegend()

p1 <- VlnPlot(sobj.unperturbed, features = c("Cd4")) + ggtitle("Unperturbed - Cd4")
p2 <- VlnPlot(sobj.tox2, features = c("Cd4")) + ggtitle("Tox2 - Cd4")
p3 <- VlnPlot(sobj.arid5b, features = c("Cd4")) + ggtitle("Arid5b - Cd4")
p4 <- VlnPlot(sobj.dvl2, features = c("Cd4")) + ggtitle("Dvl2 - Cd4")
(p1|p2|p3|p4) & NoLegend()

p1 <- VlnPlot(sobj.unperturbed, features = c("Cd8a")) + ggtitle("Unperturbed - Cd8a")
p2 <- VlnPlot(sobj.tox2, features = c("Cd8a")) + ggtitle("Tox2 - Cd8a")
p3 <- VlnPlot(sobj.arid5b, features = c("Cd8a")) + ggtitle("Arid5b - Cd8a")
p4 <- VlnPlot(sobj.dvl2, features = c("Cd8a")) + ggtitle("Dvl2 - Cd8a")
(p1|p2|p3|p4) & NoLegend()

p1 <- VlnPlot(sobj.unperturbed, features = c("Lck")) + ggtitle("Unperturbed - Lck")
p2 <- VlnPlot(sobj.tox2, features = c("Lck")) + ggtitle("Tox2 - Lck")
p3 <- VlnPlot(sobj.arid5b, features = c("Lck")) + ggtitle("Arid5b - Lck")
p4 <- VlnPlot(sobj.dvl2, features = c("Lck")) + ggtitle("Dvl2 - Lck")
(p1|p2|p3|p4) & NoLegend()

FeaturePlot(sobj.unperturbed, features = c("Cd3e","Cd4","Cd8a","Lck"), ncol = 4) & NoLegend()
FeaturePlot(sobj.tox2, features = c("Cd3e","Cd4","Cd8a","Lck"), ncol = 4) & NoLegend()
FeaturePlot(sobj.arid5b, features = c("Cd3e","Cd4","Cd8a","Lck"), ncol = 4) & NoLegend()
FeaturePlot(sobj.dvl2, features = c("Cd3e","Cd4","Cd8a","Lck"), ncol = 4) & NoLegend()

# There are several marker genes that can distinguish different types of T-cells in mice.
# Cd4, Cd8a, Foxp3, Gata3, Tbet, Rorc
p1 <- VlnPlot(sobj.unperturbed, features = c("Cd4")) + ggtitle("Unperturbed - Cd4")
p2 <- VlnPlot(sobj.tox2, features = c("Cd4")) + ggtitle("Tox2 - Cd4")
p3 <- VlnPlot(sobj.arid5b, features = c("Cd4")) + ggtitle("Arid5b - Cd4")
p4 <- VlnPlot(sobj.dvl2, features = c("Cd4")) + ggtitle("Dvl2 - Cd4")
(p1|p2|p3|p4) & NoLegend()

p1 <- VlnPlot(sobj.unperturbed, features = c("Cd8a")) + ggtitle("Unperturbed - Cd8a")
p2 <- VlnPlot(sobj.tox2, features = c("Cd8a")) + ggtitle("Tox2 - Cd8a")
p3 <- VlnPlot(sobj.arid5b, features = c("Cd8a")) + ggtitle("Arid5b - Cd8a")
p4 <- VlnPlot(sobj.dvl2, features = c("Cd8a")) + ggtitle("Dvl2 - Cd8a")
(p1|p2|p3|p4) & NoLegend()

p1 <- VlnPlot(sobj.unperturbed, features = c("Foxp3")) + ggtitle("Unperturbed - Foxp3")
p2 <- VlnPlot(sobj.tox2, features = c("Foxp3")) + ggtitle("Tox2 - Foxp3")
p3 <- VlnPlot(sobj.arid5b, features = c("Foxp3")) + ggtitle("Arid5b - Foxp3")
p4 <- VlnPlot(sobj.dvl2, features = c("Foxp3")) + ggtitle("Dvl2 - Foxp3")
(p1|p2|p3|p4) & NoLegend()

p1 <- VlnPlot(sobj.unperturbed, features = c("Gata3")) + ggtitle("Unperturbed - Gata3")
p2 <- VlnPlot(sobj.tox2, features = c("Gata3")) + ggtitle("Tox2 - Gata3")
p3 <- VlnPlot(sobj.arid5b, features = c("Gata3")) + ggtitle("Arid5b - Gata3")
p4 <- VlnPlot(sobj.dvl2, features = c("Gata3")) + ggtitle("Dvl2 - Gata3")
(p1|p2|p3|p4) & NoLegend()

p1 <- VlnPlot(sobj.unperturbed, features = c("Tbet")) + ggtitle("Unperturbed - Tbet")
p2 <- VlnPlot(sobj.tox2, features = c("Tbet")) + ggtitle("Tox2 - Tbet")
p3 <- VlnPlot(sobj.arid5b, features = c("Tbet")) + ggtitle("Arid5b - Tbet")
p4 <- VlnPlot(sobj.dvl2, features = c("Tbet")) + ggtitle("Dvl2 - Tbet")
(p1|p2|p3|p4) & NoLegend()

p1 <- VlnPlot(sobj.unperturbed, features = c("Rorc")) + ggtitle("Unperturbed - Rorc")
p2 <- VlnPlot(sobj.tox2, features = c("Rorc")) + ggtitle("Tox2 - Rorc")
p3 <- VlnPlot(sobj.arid5b, features = c("Rorc")) + ggtitle("Arid5b - Rorc")
p4 <- VlnPlot(sobj.dvl2, features = c("Rorc")) + ggtitle("Dvl2 - Rorc")
(p1|p2|p3|p4) & NoLegend()

FeaturePlot(sobj.unperturbed, features = c("Cd4","Cd8a","Foxp3","Gata3"), ncol = 4) & NoLegend()
FeaturePlot(sobj.tox2, features = c("Cd4","Cd8a","Foxp3","Gata3"), ncol = 4) & NoLegend()
FeaturePlot(sobj.arid5b, features = c("Cd4","Cd8a","Foxp3","Gata3"), ncol = 4) & NoLegend()
FeaturePlot(sobj.dvl2, features = c("Cd4","Cd8a","Foxp3","Gata3"), ncol = 4) & NoLegend()

# Knockout 샘플 별 전체 세포수 대비 state 별 세포수 비율 비교 ----
as_tibble(metadata) %>% 
  group_by(condition, state) %>% 
  summarise(cell_n = n()) %>%
  add_count(condition, wt = cell_n, name = "total_cell") %>% 
  mutate(cell_p = cell_n/total_cell) %>% 
  arrange(desc(total_cell)) -> metadata.summary

p0 <- metadata.summary %>% 
  ggscatter(x = "total_cell", y = "cell_p", add = "reg.line", 
            add.params = list(color = "blue", fill = "lightgray"), conf.int = T) +
  stat_cor(method = "pearson", label.x = 3000) + ggtitle("all")
p1 <- metadata.summary %>% 
  filter(state == "cycling") %>% 
  ggscatter(x = "total_cell", y = "cell_p", add = "reg.line", 
            add.params = list(color = "blue", fill = "lightgray"), conf.int = T) +
  stat_cor(method = "pearson", label.x = 3000, label.y = 0.7) + ggtitle("cycling")
p2 <- metadata.summary %>% 
  filter(state == "effector") %>% 
  ggscatter(x = "total_cell", y = "cell_p", add = "reg.line", 
            add.params = list(color = "blue", fill = "lightgray"), conf.int = T) +
  stat_cor(method = "pearson", label.x = 3000) + ggtitle("effector")
p3 <- metadata.summary %>% 
  filter(state == "other") %>% 
  ggscatter(x = "total_cell", y = "cell_p", add = "reg.line", 
            add.params = list(color = "blue", fill = "lightgray"), conf.int = T) +
  stat_cor(method = "pearson", label.x = 3000) + ggtitle("other")
p4 <- metadata.summary %>% 
  filter(state == "progenitor") %>% 
  ggscatter(x = "total_cell", y = "cell_p", add = "reg.line", 
            add.params = list(color = "blue", fill = "lightgray"), conf.int = T) +
  stat_cor(method = "pearson", label.x = 3000) + ggtitle("progenitor")
p5 <- metadata.summary %>% 
  filter(state == "terminal exhausted") %>% 
  ggscatter(x = "total_cell", y = "cell_p", add = "reg.line", 
            add.params = list(color = "blue", fill = "lightgray"), conf.int = T) +
  stat_cor(method = "pearson", label.x = 3000, label.y = 0.6) + ggtitle("terminal exhausted")
ggarrange(p0,p1,p2,p3,p4,p5, ncol = 3, nrow = 2)

# 'other' state에 대한 조사 ----
sobj.markers <- FindAllMarkers(sobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sobj.markers %>% 
  filter(cluster == 11) %>% 
  slice_max(n = 10, order_by = avg_log2FC) %>% 
  select(gene) %>% as_vector() -> markers
sobj.markers %>% 
  group_by(cluster) %>% 
  slice_max(n = 10, order_by = avg_log2FC) %>% 
  ungroup() %>% select(gene) %>% as_vector() -> markers
VlnPlot(sobj, features = markers, group.by = "state", ncol = 5)
FeaturePlot(sobj, features = markers, ncol = 5)
sobj.markers %>% 
  group_by(cluster) %>% 
  top_n(n=10, wt=avg_log2FC) -> top10
DoHeatmap(sobj, features = top10$gene, group.by = "state", size = 3)

mat <- sobj[["RNA"]]@data[top10$gene,] %>% as.matrix()
mat <- t(scale(t(mat)))
cluster_anno <- sobj@meta.data$state
quantile(mat, c(0.1, 0.95))
col_fun = circlize::colorRamp2(c(-1, 0, 3), c("#FF00FF", "black", "#FFFF00"))
Heatmap(mat, name = "Expression   ",  
        column_split = factor(cluster_anno),
        cluster_columns = TRUE,
        show_column_dend = F,
        cluster_column_slices = TRUE,
        column_title_gp = gpar(fontsize = 10),
        column_gap = unit(0.5, "mm"),
        cluster_rows = TRUE,
        show_row_dend = F,
        col = col_fun,
        row_names_gp = gpar(fontsize = 10),
        column_title_rot = 45,
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))),
        show_column_names = FALSE,
        use_raster = TRUE,
        raster_quality = 4)
# Unpurterbed에서 66 genes의 stacked violin plot 생성 ----
# The following requested variables were not found: Fzd1, P2rx7 -> 2 genes removed!
# Warning: All cells have the same value of Fzd3 -> Removed!
as_tibble(sobj@meta.data) %>% 
  filter(condition != "Unperturbed", condition != "Fzd1", condition != "P2rx7", condition != "Fzd3") %>% 
  distinct(condition) %>% as_vector() -> markers
p1 <- VlnPlot(sobj.unperturbed, features = markers[c(1:21)], stack = T, flip = T, assay = "RNA") + NoLegend() + ggtitle("1-21")
p2 <- VlnPlot(sobj.unperturbed, features = markers[c(22:42)], stack = T, flip = T, assay = "RNA") + NoLegend() + ggtitle("22-42")
p3 <- VlnPlot(sobj.unperturbed, features = markers[c(43:63)], stack = T, flip = T, assay = "RNA") + NoLegend() + ggtitle("43-63")
p1|p2|p3

as_tibble(sobj@meta.data) %>% 
  filter(condition != "Unperturbed", condition != "Fzd1", condition != "P2rx7", condition != "Fzd3") %>% 
  distinct(condition) %>% as_vector() -> markers
p1 <- VlnPlot(sobj.unperturbed, features = markers[c(1:21)], stack = T, flip = T, assay = "RNA", group.by = "state") + NoLegend() + ggtitle("1-21")
p2 <- VlnPlot(sobj.unperturbed, features = markers[c(1:21)], stack = T, flip = T, assay = "RNA", group.by = "state") + NoLegend() + ggtitle("22-42")
p3 <- VlnPlot(sobj.unperturbed, features = markers[c(1:21)], stack = T, flip = T, assay = "RNA", group.by = "state") + NoLegend() + ggtitle("43-63")
p1|p2|p3

RidgePlot(sobj.unperturbed, features = markers[c(1:21)], ncol = 3, group.by = "state")
DotPlot(sobj.unperturbed, features = markers[c(1:21)], group.by = "state")
DoHeatmap(sobj.unperturbed, features = markers, group.by = "state")

# ComplexHeatmap
mat <- sobj.unperturbed[["RNA"]]@data[markers,] %>% as.matrix()
mat <- t(scale(t(mat)))
cluster_anno <- sobj.unperturbed@meta.data$state
quantile(mat, c(0.1, 0.95))
col_fun = circlize::colorRamp2(c(-1, 0, 3), c("#FF00FF", "black", "#FFFF00"))
Heatmap(mat, name = "Expression   ",  
        column_split = factor(cluster_anno),
        cluster_columns = TRUE,
        show_column_dend = F,
        cluster_column_slices = TRUE,
        column_title_gp = gpar(fontsize = 10),
        column_gap = unit(0.5, "mm"),
        cluster_rows = TRUE,
        show_row_dend = F,
        col = col_fun,
        row_names_gp = gpar(fontsize = 10),
        column_title_rot = 45,
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))),
        show_column_names = FALSE,
        use_raster = TRUE,
        raster_quality = 4)

# progenitor와 terminal exhausted 간 발현량 비교 통계 검정
metadata %>% 
  count(condition) %>% arrange(desc(n)) %>% 
  filter(n<100, condition != "Unperturbed", condition != "Fzd1", condition != "P2rx7", condition != "Fzd3") %>% 
  select(condition) %>% as_vector() -> markers.lowcell
my_comparisons <- list(c("progenitor", "terminal exhausted"))
VlnPlot(sobj.unperturbed, features = markers.lowcell[c(1:10)], group.by = "state", ncol = 5)
VlnPlot(sobj.unperturbed, features = markers.lowcell[1], group.by = "state") + NoLegend() +
  stat_compare_means(comparisons = my_comparisons, label.y = 3) + ylim(0, 4)
for(i in 21:29) {
  assign(paste0("p",i), 
         VlnPlot(sobj.unperturbed, features = markers.lowcell[i], group.by = "state") + NoLegend() +
           stat_compare_means(comparisons = my_comparisons, label.y = 3) + ylim(0, 4)
  )
}
p21+p22|p23+p24|p25+p26|p27+p28|p29+p29

# More common expressed genes ----
p1 <- VlnPlot(sobj.unperturbed, features = c("Cd8b1"))
p2 <- VlnPlot(sobj.tox2, features = c("Cd8b1"))
p3 <- VlnPlot(sobj.arid5b, features = c("Cd8b1"))
p4 <- VlnPlot(sobj.dvl2, features = c("Cd8b1"))
(p1|p2|p3|p4) & NoLegend()

# ProjectTILs ----
mouse.cd8 <- read_rds("/db2/users/shnam/04.IMP/02.data/ProjectTILs_reference_maps/ref_LCMV_Atlas_mouse_v1.rds")
mouse.cd4 <- read_rds("/db2/users/shnam/04.IMP/02.data/ProjectTILs_reference_maps/ref_LCMV_CD4_mouse_release_v1.rds")
mouse.til <- read_rds("/db2/users/shnam/04.IMP/02.data/ProjectTILs_reference_maps/ref_TILAtlas_mouse_v1.rds")

refCols <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000", "#87f6a5", "#e812dd")
p1 <- DimPlot(mouse.til, reduction = "umap", label = T, cols = refCols) + ggtitle("Reference map") + NoLegend()

markers <- c("Cd4", "Cd8a", "Ccr7", "Tcf7", "Pdcd1", "Havcr2", "Tox", "Izumo1r", "Cxcr6", "Xcl1", "Gzmb", "Gzmk", "Ifng", "Foxp3")
VlnPlot(mouse.til, features = markers, stack = T, flip = T, assay = "RNA") + NoLegend() + ggtitle("Reference")
p3 <- VlnPlot(sobj, features = markers, stack = T, flip = T, assay = "RNA") + NoLegend()

query <- Run.ProjecTILs(sobj, ref = mouse.til)
p2 <- plot.projection(mouse.til, query, linesize = 0.5, pointsize = 0.5) + NoLegend()
p1|p2
p4 <- plot.statepred.composition(mouse.til, query, metric = "Percent")
p3|p4

plot.states.radar(mouse.til, query = query, min.cells = 30)

p1 <- DimPlot(sobj, group.by = "state", label = T) + NoLegend()
querydata <- ProjecTILs.classifier(query = sobj, ref = mouse.til)
palette <- c("#FF0000", "#00B6EB", "#53B400", "#F8766D", "#A58AFF", "#d1cfcc", "#edbe2a",
             "#87f6a5", "#e812dd", "#777777")
names(palette) <- c(levels(mouse.til$functional.cluster), "NA")
p2 <- DimPlot(querydata, group.by = "functional.cluster", cols = palette)
p1|p2

# 대희님 요청사항----
# Unperturbed와 Pdcd1/Tox knockout 샘플 간 T cell state, cell number 비교를 통해 exhausted 세포수가 유의하게 적은지 확인
metadata.summary %>% 
  filter(condition == "Unperturbed" | condition == "Tox" | condition == "Pdcd1")
metadata.summary %>% 
  filter(state == "terminal exhausted") %>% 
  ggscatter(x = "total_cell", y = "cell_n", add = "reg.line", label = "condition", repel = T,
            label.select = c("Unperturbed", "Tox"),
            add.params = list(color = "blue", fill = "lightgray"), conf.int = T) +
  stat_cor(method = "pearson", label.x = 3000, label.y = 0.6) + ggtitle("terminal exhausted") -> p1
metadata.summary %>% 
  filter(state == "terminal exhausted") %>% 
  ggscatter(x = "total_cell", y = "cell_p", add = "reg.line", label = "condition", repel = T,
            label.select = c("Unperturbed", "Tox"),
            add.params = list(color = "blue", fill = "lightgray"), conf.int = T) +
  stat_cor(method = "pearson", label.x = 3000, label.y = 0.6) + ggtitle("terminal exhausted") -> p2
p1|p2


