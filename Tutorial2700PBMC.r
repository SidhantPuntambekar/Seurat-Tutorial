library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "~/seuratTutorial/data/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
head(pbmc@meta.data, 5)

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Perform linear dimension reduction
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

#Determine Dimensionality of data set. May take time to run so comment out if taking too much CPU
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)

#Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10) #Construct KNN graph and refine edge weights between two cells based on the shared overlap of local neighbors using Jaccard Index
pbmc <- FindClusters(pbmc, resolution = 0.5) #Run Louvain algorithm to cluster cells together with resolution parameter of 0.5. Useful for good returns of single cell datasets around 3000 cells
head(Idents(pbmc), 5) #Look at cluster IDs of the first five cells.

#Run non-linear dimensional reduction (UMAP/tSNE)
pbmc <- RunUMAP(pbmc, dims = 1:10) #RunUMAP vis tool on pbmc with dims of 1:10
DimPlot(pbmc, reduction = "umap") #Create DimPlot to visualize individual clusters

#Finding differentially expressed features (cluster biomarkers)
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25) #Find all markers for differential expression using a min percentage in either two groups of cells
head(cluster1.markers, n = 5) #Check first five markers of cluster 1
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25) #Find all markers distinguishing cluster 5 from clusters 0 and 3
head(cluster5.markers, n = 5) #Check first five markers of cluster 5
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #Find markers for every cluster compared to all remaining cells report only the positve ones
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE) #Conduct ROC test on cluster 1. Returns classification power for any individual marker
VlnPlot(pbmc, features = c("MS4A1", "CD79A")) #Construct violin plot to show expression probability distributions across clusters
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE) #Plot raw counts as well for two different cells
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A")) #Create feature plot of nine cell types for UMAPs
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) #Create top 10 markers for each cluster
DoHeatmap(pbmc, features = top10$gene) + NoLegend() #Create heat map of top 10 markers

#Assigning cell type identity to clusters
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet") #Assign cell types to clusters
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
head(Idents(pbmc), 5) # Look at cluster IDs of the first 5 cells