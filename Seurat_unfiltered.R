# Load necessary libraries
#.libPaths(c("/cluster/home/petear/R", .libPaths())) #commented when running locally
library("dplyr")
library("Seurat")
library("patchwork")
library("ggplot2")

# Create directories for unfiltered analysis
dir.create("unfiltered_analysis", showWarnings = FALSE)
dir.create("unfiltered_figures", showWarnings = FALSE)

# Load data
#pbmc.data <- Read10X(data.dir = "/cluster/projects/nn4605k/peter/from_junbai/to_peter_data/scRNAseq/NT8")
pbmc.data <- Read10X(data.dir = "/home/petear/sc_seq_data/to_peter_data/scRNAseq/NT8")

# Create UNFILTERED Seurat object
pbmc.unfiltered <- CreateSeuratObject(counts = pbmc.data, project = "JBW_master_unfiltered",
                                     min.cells = 0, min.features = 0)

# Calculate mitochondrial percentage
pbmc.unfiltered[["percent.mt"]] <- PercentageFeatureSet(pbmc.unfiltered, pattern = "^MT-")

# Create comprehensive PDF for UNFILTERED analysis
generate_plots <- function(pdf_file, width, height) {
  pdf(pdf_file, width = width, height = height)

# Sink output to a text file for logging
sink("unfiltered_analysis/unfiltered_analysis_log.txt")

# 1. QUALITY CONTROL METRICS - UNFILTERED
cat("=== UNFILTERED DATA SUMMARY ===\n")
cat("Total cells:", ncol(pbmc.unfiltered), "\n")
cat("Total genes:", nrow(pbmc.unfiltered), "\n")
cat("Range of genes per cell:", min(pbmc.unfiltered$nFeature_RNA), "-", max(pbmc.unfiltered$nFeature_RNA), "\n")
cat("Range of UMIs per cell:", min(pbmc.unfiltered$nCount_RNA), "-", max(pbmc.unfiltered$nCount_RNA), "\n")
cat("Range of mitochondrial %:", min(pbmc.unfiltered$percent.mt), "-", max(pbmc.unfiltered$percent.mt), "\n")

# Violin plots of QC metrics
vln_plot <- VlnPlot(pbmc.unfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                   ncol = 3, pt.size = 0.1) +
  ggtitle("QC Violin Plots - UNFILTERED DATA")
print(vln_plot)

# Scatter plots showing relationships
plot1 <- FeatureScatter(pbmc.unfiltered, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  ggtitle("UMIs vs Mitochondrial % - UNFILTERED") +
  theme(legend.position = "none")
plot2 <- FeatureScatter(pbmc.unfiltered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  ggtitle("UMIs vs Genes - UNFILTERED") +
  theme(legend.position = "none")
combined_scatter <- plot1 + plot2 + plot_annotation(title = "QC Scatter Plots - UNFILTERED DATA")
print(combined_scatter)

# Distribution histograms
dist1 <- ggplot(pbmc.unfiltered@meta.data, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  labs(title = "Distribution of Genes per Cell - UNFILTERED",
       subtitle = paste("Median:", median(pbmc.unfiltered$nFeature_RNA), "genes"),
       x = "Number of Genes", y = "Number of Cells") +
  theme_minimal()
dist2 <- ggplot(pbmc.unfiltered@meta.data, aes(x = nCount_RNA)) +
  geom_histogram(bins = 50, fill = "forestgreen", alpha = 0.7) +
  labs(title = "Distribution of UMIs per Cell - UNFILTERED",
       subtitle = paste("Median:", median(pbmc.unfiltered$nCount_RNA), "UMIs"),
       x = "Number of UMIs", y = "Number of Cells") +
  theme_minimal()
dist3 <- ggplot(pbmc.unfiltered@meta.data, aes(x = percent.mt)) +
  geom_histogram(bins = 50, fill = "coral", alpha = 0.7) +
  labs(title = "Distribution of Mitochondrial % - UNFILTERED",
       subtitle = paste("Median:", round(median(pbmc.unfiltered$percent.mt), 2), "%"),
       x = "Mitochondrial %", y = "Number of Cells") +
  theme_minimal()
print(dist1)
print(dist2)
print(dist3)

# 2. NORMALIZATION AND VARIABLE FEATURES - UNFILTERED
cat("\n=== NORMALIZING UNFILTERED DATA ===\n")
pbmc.unfiltered <- NormalizeData(pbmc.unfiltered, normalization.method = "LogNormalize",
                                scale.factor = 10000)
pbmc.unfiltered <- FindVariableFeatures(pbmc.unfiltered, selection.method = "vst",
                                       nfeatures = 2000)

# Identify top variable genes
top10_unfiltered <- head(VariableFeatures(pbmc.unfiltered), 10)
cat("Top 10 variable genes in unfiltered data:\n")
print(top10_unfiltered)

# Plot variable features
var_plot1 <- VariableFeaturePlot(pbmc.unfiltered) +
  ggtitle("Variable Features - UNFILTERED DATA")
var_plot2 <- LabelPoints(plot = var_plot1, points = top10_unfiltered, repel = TRUE) +
  ggtitle("Top 10 Variable Genes - UNFILTERED DATA")
print(var_plot1)
print(var_plot2)

# 3. DATA SCALING - UNFILTERED
cat("\n=== SCALING UNFILTERED DATA ===\n")
all.genes.unfiltered <- rownames(pbmc.unfiltered)
pbmc.unfiltered <- ScaleData(pbmc.unfiltered, features = all.genes.unfiltered)

# 4. PCA - UNFILTERED
cat("\n=== PERFORMING PCA ON UNFILTERED DATA ===\n")
pbmc.unfiltered <- RunPCA(pbmc.unfiltered, features = VariableFeatures(object = pbmc.unfiltered))

# PCA visualizations
print(pbmc.unfiltered[["pca"]], dims = 1:5, nfeatures = 5)
pca_loadings <- VizDimLoadings(pbmc.unfiltered, dims = 1:2, reduction = "pca") +
  ggtitle("PCA Loadings (PC1 & PC2) - UNFILTERED DATA")
print(pca_loadings)
pca_plot <- DimPlot(pbmc.unfiltered, reduction = "pca") +
  NoLegend() +
  ggtitle("PCA Visualization - UNFILTERED DATA")
print(pca_plot)

# PCA heatmaps with proper titles
print_dimheatmap_with_title <- function(seurat_obj, dims, title) {
  heatmap_plot <- DimHeatmap(seurat_obj, dims = dims, cells = 500, balanced = TRUE, fast = FALSE)
  print(heatmap_plot)
  grid::grid.text(title, x = 0.5, y = 0.95, gp = grid::gpar(fontsize = 16, fontface = "bold"))
}

print_dimheatmap_with_title(pbmc.unfiltered, dims = 1, "PCA Heatmap - PC1 - UNFILTERED DATA")
print_dimheatmap_with_title(pbmc.unfiltered, dims = 1:15, "PCA Heatmap - PCs 1-15 - UNFILTERED DATA")

# Elbow plot
elbow_plot <- ElbowPlot(pbmc.unfiltered) +
  ggtitle("PCA Elbow Plot - UNFILTERED DATA")
print(elbow_plot)

# 5. CLUSTERING - UNFILTERED
cat("\n=== CLUSTERING UNFILTERED DATA ===\n")
pbmc.unfiltered <- FindNeighbors(pbmc.unfiltered, dims = 1:10)
pbmc.unfiltered <- FindClusters(pbmc.unfiltered, resolution = 0.5)

# Check cluster sizes
cluster_sizes <- table(Idents(pbmc.unfiltered))
cat("Cluster sizes in unfiltered data:\n")
print(cluster_sizes)

# 6. UMAP - UNFILTERED
cat("\n=== RUNNING UMAP ON UNFILTERED DATA ===\n")
pbmc.unfiltered <- RunUMAP(pbmc.unfiltered, dims = 1:10)

# UMAP visualization
umap_plot <- DimPlot(pbmc.unfiltered, reduction = "umap", label = TRUE) +
  ggtitle("UMAP Clustering - UNFILTERED DATA")
print(umap_plot)

# 7. DIFFERENTIAL EXPRESSION - UNFILTERED
cat("\n=== FINDING MARKERS IN UNFILTERED DATA ===\n")
pbmc.markers.unfiltered <- FindAllMarkers(pbmc.unfiltered, only.pos = TRUE,
                                         min.pct = 0.25, logfc.threshold = 0.25)

# Summary of markers found
marker_summary <- pbmc.markers.unfiltered %>%
  group_by(cluster) %>%
  summarise(n_markers = n(),
            top_marker = first(gene),
            avg_log2FC = round(first(avg_log2FC), 2))
cat("Marker summary by cluster:\n")
print(marker_summary)

# 8. MARKER VISUALIZATION - UNFILTERED
# Common cell type markers
common_markers <- c("MS4A1", "CD79A", "CD79B", "CD19", "CD3E", "CD3D", "CD4", "CD8A",
                   "NKG7", "GNLY", "CD14", "LYZ", "FCGR3A", "FCER1A", "PPBP")

# Check which markers are present in our dataset
available_markers <- common_markers[common_markers %in% rownames(pbmc.unfiltered)]
cat("Available common markers in dataset:\n")
print(available_markers)

# Feature plots for available markers
if(length(available_markers) > 0) {
  feature_plot <- FeaturePlot(pbmc.unfiltered, features = available_markers[1:min(9, length(available_markers))]) +
    ggtitle("Common Marker Expression - UNFILTERED DATA")
  print(feature_plot)
}

# Violin plots for key markers
if("MS4A1" %in% rownames(pbmc.unfiltered) & "CD79A" %in% rownames(pbmc.unfiltered)) {
  vln_markers <- VlnPlot(pbmc.unfiltered, features = c("MS4A1", "CD79A")) +
    ggtitle("B Cell Markers - UNFILTERED DATA")
  print(vln_markers)
}

# Heatmap of top markers
top_markers_unfiltered <- pbmc.markers.unfiltered %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
heatmap_plot <- DoHeatmap(pbmc.unfiltered, features = top_markers_unfiltered$gene) +
  ggtitle("Top 5 Markers per Cluster - UNFILTERED DATA")
print(heatmap_plot)

# 9. QUALITY ASSESSMENT BY CLUSTER - UNFILTERED
# Add QC metrics to UMAP
feature_umi <- FeaturePlot(pbmc.unfiltered, features = "nCount_RNA") +
  ggtitle("UMI Count per Cell - UNFILTERED DATA")
feature_genes <- FeaturePlot(pbmc.unfiltered, features = "nFeature_RNA") +
  ggtitle("Gene Count per Cell - UNFILTERED DATA")
feature_mt <- FeaturePlot(pbmc.unfiltered, features = "percent.mt") +
  ggtitle("Mitochondrial % per Cell - UNFILTERED DATA")
print(feature_umi)
print(feature_genes)
print(feature_mt)

# Violin plots of QC metrics by cluster
qc_by_cluster <- VlnPlot(pbmc.unfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                        group.by = "seurat_clusters", ncol = 3) +
  ggtitle("QC Metrics by Cluster - UNFILTERED DATA")
print(qc_by_cluster)

###Changed
# Filtering step
pbmc.filtered <- subset(pbmc.unfiltered, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10)

# Additional violin plots for more details
vln_plot_genes <- VlnPlot(pbmc.unfiltered, features = c("nFeature_RNA"), ncol = 1, pt.size = 0.1) +
  ggtitle("Detailed Violin Plot for Genes per Cell")
print(vln_plot_genes)

vln_plot_umis <- VlnPlot(pbmc.unfiltered, features = c("nCount_RNA"), ncol = 1, pt.size = 0.1) +
  ggtitle("Detailed Violin Plot for UMIs per Cell")
print(vln_plot_umis)

vln_plot_mt <- VlnPlot(pbmc.unfiltered, features = c("percent.mt"), ncol = 1, pt.size = 0.1) +
  ggtitle("Detailed Violin Plot for Mitochondrial %")
print(vln_plot_mt)

# Close PDF device and stop sinking output
dev.off()
sink()
}

generate_plots("unfiltered_figures/RNA_seq_UNFILTERED_analysis_BIG.pdf", 25, 10)
generate_plots("unfiltered_figures/RNA_seq_UNFILTERED_analysis.pdf", 15, 10)

# 10. EXPORT UNFILTERED DATA AND RESULTS
# Save Seurat object
saveRDS(pbmc.unfiltered, "unfiltered_analysis/seurat_object_unfiltered.rds")

# Export metadata
write.csv(pbmc.unfiltered@meta.data, "unfiltered_analysis/metadata_unfiltered.csv")

# Export differential expression results
write.csv(pbmc.markers.unfiltered, "unfiltered_analysis/differential_expression_unfiltered.csv")

# Export cluster markers summary
write.csv(marker_summary, "unfiltered_analysis/cluster_markers_summary_unfiltered.csv")

# Export variable features
variable_features <- VariableFeatures(pbmc.unfiltered)
write.csv(data.frame(Gene = variable_features), "unfiltered_analysis/variable_features_unfiltered.csv")

# 11. CREATE SUMMARY REPORT
sink("unfiltered_analysis/unfiltered_analysis_summary.txt")
cat("UNFILTERED SINGLE-CELL RNA-SEQ ANALYSIS SUMMARY\n")
cat("===============================================\n\n")
cat("Dataset: JB_master_unfiltered\n")
cat("Analysis date:", date(), "\n\n")
cat("DATA OVERVIEW:\n")
cat("-------------\n")
cat("Total cells:", ncol(pbmc.unfiltered), "\n")
cat("Total genes:", nrow(pbmc.unfiltered), "\n")
cat("Median genes per cell:", median(pbmc.unfiltered$nFeature_RNA), "\n")
cat("Median UMIs per cell:", median(pbmc.unfiltered$nCount_RNA), "\n")
cat("Median mitochondrial %:", round(median(pbmc.unfiltered$percent.mt), 2), "%\n\n")
cat("CLUSTERING RESULTS:\n")
cat("------------------\n")
cat("Number of clusters:", length(unique(Idents(pbmc.unfiltered))), "\n")
cat("Cluster sizes:\n")
print(cluster_sizes)
cat("\n")
cat("VARIABLE FEATURES:\n")
cat("------------------\n")
cat("Number of variable features:", length(variable_features), "\n")
cat("Top 10 variable genes:\n")
print(top10_unfiltered)
cat("\n")
cat("MARKER DISCOVERY:\n")
cat("----------------\n")
cat("Total markers found:", nrow(pbmc.markers.unfiltered), "\n")
cat("Markers per cluster:\n")
print(marker_summary)
sink()

# Final message
cat("\n=== UNFILTERED ANALYSIS COMPLETE ===\n")
cat("Files saved in:\n")
cat("- unfiltered_analysis/: All data files and results\n")
cat("- unfiltered_figures/: All visualization PDFs\n")
cat("\nKey files:\n")
cat("- seurat_object_unfiltered.rds: Complete Seurat object\n")
cat("- metadata_unfiltered.csv: Cell metadata with clusters\n")
cat("- differential_expression_unfiltered.csv: All marker genes\n")
cat("- RNA_seq_UNFILTERED_analysis.pdf: All visualizations\n")
cat("- unfiltered_analysis_summary.txt: Analysis summary\n")

# Print final summary to console
cat("\nFinal unfiltered dataset statistics:\n")
cat("Cells:", ncol(pbmc.unfiltered), "\n")
cat("Genes:", nrow(pbmc.unfiltered), "\n")
cat("Clusters:", length(unique(Idents(pbmc.unfiltered))), "\n")
cat("Analysis ready for comparison with filtered dataset!\n")
