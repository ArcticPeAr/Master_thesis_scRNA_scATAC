library(Seurat)
library(ggplot2)
library(patchwork)

# Get input/output from Snakemake
input_rds <- snakemake@input[["rds"]]
output_png <- snakemake@output[["qc_plot"]]
sample <- snakemake@wildcards$sample
version <- snakemake@wildcards$version  # "unfiltered" or "filtered"

# Load Seurat object
pbmc <- readRDS(input_rds)

# Create QC plots (simplified version of your original)
pdf_file <- paste0("figures/", sample, "/", version, "_qc_summary.pdf")
pdf(pdf_file, width = 15, height = 10)

# Violin plots
vln_plot <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                   ncol = 3, pt.size = 0.1) +
  ggtitle(paste(sample, version, "- QC Violin Plots"))
print(vln_plot)

# Scatter plots
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
combined_scatter <- plot1 + plot2 + 
  plot_annotation(title = paste(sample, version, "- QC Scatter Plots"))
print(combined_scatter)

# Histograms
dist1 <- ggplot(pbmc@meta.data, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 50, fill = "steelblue") +
  labs(title = "Genes per Cell Distribution")
print(dist1)

dev.off()

# Save single PNG for Snakemake output
ggsave(output_png, vln_plot, width = 12, height = 4)