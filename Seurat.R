#Script based on Seurat v5 tutorial at https://satijalab.org/seurat/articles/pbmc3k_tutorial.html and https://satijalab.org/seurat/articles/visualization_vignette


# Sequencing: to get raw reads of the gene sequences and cellular barcodes ( and potentially UMI ).
# Demultiplexing (e.g., using Cell Ranger): cellular barcodes to assign reads to their cell of origin.
# Alignment: Map the gene sequences from the reads to a reference genome/transcriptome.
# Quantification: Count the number of molecules (e.g., UMIs) detected for each gene in each cell.
# Creates the count matrix (.mtx.gz) and associated files listing genes (features.tsv.gz - also need the refeerence for this to name genes) and cell barcodes (.tsv.gz).




.libPaths(c("/cluster/home/user/R", .libPaths())) 
library("dplyr")
library("Seurat")
library("patchwork")
library("ggplot2")  

#? RNA seq 
## Load data in the form of 
# barcodes.tsv.gz
# features.tsv.gz
# matrix.mtx.gz
# pbmc_granulocyte_sorted_10k_per_barcode_metrics.csv
pbmc.data <- Read10X(data.dir = "/cluster/projects/nn4605k/peter/from_junbai/to_peter_data/scRNAseq/NT8")

#! Filtering cells and creating Seurat object
# min.cells = 3: Only include genes that are expressed in at least 3 cells
# min.features = 200: Only include cells that have at least 200 detected genes
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "JB_master", min.cells = 3, min.features = 200)


# Create PDF of QC plots
pdf("RNA_seq_QC.pdf")

#! QC here is metric on how much mitochondrial genes are found. A high number of mitochondrial genes may indicate cell damage where cytoplasmic RNA has leaked out of the cell
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#! Normalizing the data 
# Scale factor of 10,000: counts per 10,000
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Authors of Seurat recommend using SCTransform for normalization
# Normalize, scale, and find variable features using the recommended SCTransform method
# Regress out the effect of mitochondrial gene percentage (percent.mt) as a covariate
# Results are stored in the "SCT" assay
#https://satijalab.org/seurat/articles/sctransform_vignette

# FindVariableFeatures: Identifies genes with high cell-to-cell variation using the variance
# stabilizing transformation (vst - variance stabilizing transformation) method, which models the mean-variance relationship.
# These 2000 selected variable genes are used for downstream dimensionality reduction (PCA)
# and clustering, as they are likely to capture the main sources of biological heterogeneity.
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
# Highly variable genes (likely to be biologically relevant) are highlighted.
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#! Scaling the data
all.genes <- rownames(pbmc) # Get all gene names as they are rows in the pbmc
pbmc <- ScaleData(pbmc, features = all.genes) #Z-score normalization 

#! Linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


# Examine and visualize PCA results a few different ways
pbmc <- FindVariableFeatures(pbmc)

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca") # Visualize the top genes contributing to PC1 and PC2

DimPlot(pbmc, reduction = "pca") + NoLegend()   
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

#! Determining the 'dimensionality' of the data
ElbowPlot(pbmc)

#! Cluster the cells 
# Here we use the first 10 PCs, as suggested by the elbow plot
# The FindNeighbors function constructs a KNN graph based on the Euclidean distance in PCA space.
# The FindClusters function applies a graph-based clustering algorithm (Louvain algorithm by default)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

#! Non-linear dimensional reduction (UMAP/tSNE) 

pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")

#! Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)  # find genes that are expressed in cluster 2 compared to all other clusters
head(cluster2.markers, n = 5) # show the top 5 markers for cluster 2

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)


# find markers for every cluster compared to all remaining cells, report only the positive
# ones
#For selection find good markers that can be selected. 
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)  # find all positive markers
#Group by cluster identity and filter for those with avg_log2FC > 1
pbmc.markers %>%
    group_by(cluster) %>% 
    dplyr::filter(avg_log2FC > 1)

# find markers for cluster 0 with a logfc threshold of 0.25 using the ROC test and use cluster 0 of interest
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)


FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
    "CD8A"))

pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

#! Assigning cell type identity to clusters 
#Canonical markers available? 


dev.off()

#