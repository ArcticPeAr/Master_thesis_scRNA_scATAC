library("Seurat")
library("patchwork") # check what this is for



#seq_dir is the directory pointed to by snakemake
seq_dir <- snakemake@input[["seq_dir"]]

pbmc.data <- Read10X(data.dir = seq_dir)
pbmc.unfiltered <- CreateSeuratObject(counts = pbmc.data, project = "JBW_master_unfiltered",
                                     min.cells = 0, min.features = 0)


pbmc.unfiltered[["percent.mt"]] <- PercentageFeatureSet(pbmc.unfiltered, pattern = "^MT-")

# min.cells = 3: Only include genes that are expressed in at least 3 cells
# min.features = 200: Only include cells that have at least 200 detected genes
pbmc.filtered <- CreateSeuratObject(counts = pbmc.data, project = "JBW_master", min.cells = 3, min.features = 200)

pbmc.filtered[["percent.mt"]] <- PercentageFeatureSet(pbmc.unfiltered, pattern = "^MT-")

saveRDS(pbmc.unfiltered, snakemake@output[["pbmc.unfiltered"]])
saveRDS(pbmc.filtered, snakemake@output[["pbmc.filtered"]])


