#scRNA-seq QC
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)

library(scDblFinder)

#### Load scRNA-seq
###
file2project<-c('NT7'='NT1',
		'NT8'='NT2',
		'PT5'='PT3',
		'RT3'='RT1',
		'RT4'='RT2',
		'RT6'='RT3',
		'PT1'='PT1',
		'PT2'='PT2')

# Load scRNA-seq
for (file in c("PT1","PT2","RT3","RT4","NT7","NT8","PT5","RT6")) {
     project_name=file2project[file]
     #print(project_name)
     seurat_data <- Read10X(data.dir=paste0("../GSE240112/raw/scRNAseq/",file)) 
     seurat_obj <- CreateSeuratObject(counts= seurat_data, min.features=200, min.cells=3, project=project_name)
     assign(paste0(project_name,'.seu'), seurat_obj)
     rm(seurat_data)
     rm(seurat_obj)
}

#a list of original seurat obj
seurat_list <- c('PT1' = PT1.seu, 'PT2' = PT2.seu, 'PT3' = PT3.seu,
                 'RT1' = RT1.seu, 'RT2' = RT2.seu, 'RT3' = RT3.seu,
                 'NT1' = NT1.seu, 'NT2' = NT2.seu)


####functions for filtering
finddoublet <- function(seurat_obj){
  set.seed(34)
  #compute count and logcounts
  #counts<-seurat_obj@assays$RNA$counts
  #libsizes <- colSums(counts)
  #size.factors <- libsizes/mean(libsizes)
  #logcounts_sce <- log2(t(t(counts)/size.factors) + 1)
  #seurat_obj@assays$RNA$data<-logcounts_sce
  #
  sce_obj <- as.SingleCellExperiment(seurat_obj)
  sce_obj <- scDblFinder(sce_obj)
  #compute count and logcounts
  #counts <- assay(sce_obj, "counts")
  #libsizes <- colSums(counts)
  #size.factors <- libsizes/mean(libsizes)
  #logcounts_sce <- log2(t(t(counts)/size.factors) + 1)
  
  seu_new <- as.Seurat(sce_obj,counts="counts",data="counts")
  return(seu_new)
}

nMAD <- function(x,nmads=3){
  xm <- median(x)
  md <- median(abs(x-xm))
  mads <- xm+nmads*md
  return(mads)
}


#####do filter in each obj
nmad=5
count = 1
sample.nfeaure.cut <- c()
sample.ncount.cut <- c()
sample.mt.cut <- c()
seurat_list_qc <- c()
new_seurat_list <-c()
is_gene_filter=1

for (obj in seurat_list){
  tmp_name<-names(seurat_list)[count]
  obj<-seurat_list[[tmp_name]]
  print(tmp_name)
  tmp_file=paste( gsub('"','',tmp_name),".jpg",sep="")
  print(tmp_file)
  options(bitmapType="cairo")
  jpeg(filename=tmp_file)

  # add percent.mt for qc0
  obj[['percent.mt']] <- PercentageFeatureSet(obj, pattern = "^MT-")
  # add doublet info
  obj = finddoublet(obj)
  # Add number of genes per UMI for each cell to metadata
  obj[['log10GenesPerUMI']] <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)
  # Compute percent mito ratio
  obj[['mitoRatio']] <- obj$percent.mt/100

  # nMAD cut
  nfeature.upcut <- ceiling(nMAD(obj$nFeature_RNA,nmad))
  ncount.upcut <- ceiling(nMAD(obj$nCount_RNA,nmad))
  permt.upcut <- ceiling(nMAD(obj$percent.mt,nmad))
  sample.nfeaure.cut <- c(sample.nfeaure.cut, nfeature.upcut)
  sample.ncount.cut <- c(sample.ncount.cut, ncount.upcut)
  sample.mt.cut <- c(sample.mt.cut,permt.upcut)

  # perform QC
  # Filter out low quality reads using selected thresholds - these will change with experiment
  obj.filt <- subset(obj, subset = nFeature_RNA <= nfeature.upcut & nFeature_RNA > nfeature.upcut/20)
  obj.filt <- subset(obj.filt, subset = nCount_RNA <= ncount.upcut & nCount_RNA > ncount.upcut/20)
  obj.filt <- subset(obj.filt, subset = percent.mt <= min(permt.upcut,25))
  obj.filt <- subset(obj.filt, subset = scDblFinder.class %in% c('singlet'))
  #add genesPerUMI
  obj.filt <-subset(obj.filt, subset= (log10GenesPerUMI > 0.80) )
  #obj.filt@assays$RNA@counts@Dimnames[2]

  #filter specific genes or gene names
  #all_genes<-rownames(obj.filt)
  #removed_genes1<-all_genes[str_detect(all_genes, "\\.")]
  #removed_genes2<-all_genes[str_detect(all_genes, "^LINC")]
  #all_removed_genes<-c(removed_genes1, removed_genes2)

  if ( is_gene_filter ==1) {
    print("Gene filtering ... ")
    ## Output a logical vector for every gene on whether the more than zero counts per cell
    # Extract counts
    counts <- GetAssayData(object = obj.filt, layer = "counts", assay="RNA")

    # Output a logical vector for every gene on whether the more than zero counts per cell
    nonzero <- counts > 0

    ###### Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
    #keep_genes <- Matrix::rowSums(nonzero) >= 10
    #here we remove genes with both few features and not annotated properly
    keep_genes <- ( Matrix::rowSums(nonzero) >= 10  & !str_detect(rownames(nonzero),"\\.") & !str_detect(rownames(nonzero),"^LINC") )
    #LINC0

    # Only keeping those genes expressed in more than 10 cells and passed filtering
    filtered_counts <- counts[keep_genes, ]

    # Reassign to filtered Seurat object ?? error?? obj.filt is not the same as obj.filt2?? something wrong in createSeuratObject?
    obj.filt2 <- CreateSeuratObject(counts=filtered_counts, meta.data = obj.filt@meta.data)
    obj.filt2 <- as.SingleCellExperiment(obj.filt2)
    obj.filt2 <- as.Seurat(obj.filt2,counts="counts",data="counts")

  } else {
    obj.filt2<-obj.filt	  
    print("Original Kun preprocess without gene filtering ... ")
  }

  Idents(obj.filt2) <- names(seurat_list)[count]
  print(VlnPlot(obj.filt2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  seurat_list_qc <- c(seurat_list_qc,obj.filt2)
  count = count + 1
  #update obj
  new_seurat_list <- c(new_seurat_list, obj)
  dev.off()
}

names(sample.mt.cut) <- names(seurat_list)
names(sample.ncount.cut) <- names(seurat_list)
names(sample.nfeaure.cut) <- names(seurat_list)
names(seurat_list_qc) <- names(seurat_list)
names(new_seurat_list) <- names(seurat_list)
print("Done - QC ")

# Create .RData object to load at any time
save(seurat_list_qc,file="data/seurat_list_qc.RData")
save(sample.mt.cut, file="data/sample_mt_cut.RData")
save(new_seurat_list,file="data/new_seurat_list.RData")
save(sample.ncount.cut, file="data/sample_ncount_cut.RData")
save(sample.nfeaure.cut, file="data/sample_nfeature_cut.RData")
#load(data/seurat_list_qc.RData")
print("Export data done ")

#### Create a merged Seurat object
#seurat_list <- seurat_list_qc
#merged_seurat<- merge(x=seurat_list_qc[[1]], y=seurat_list_qc[2:length(seurat_list_qc)], add.cell.id=c("PT1","PT2","PT3","RT1","RT2","RT3","NT1","NT2"))
merged_seurat<- merge(x=seurat_list_qc[[1]], y=seurat_list_qc[2:length(seurat_list_qc)], add.cell.id=names(seurat_list_qc))


# Check that the merged object has the appropriate sample-specific prefixes
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)

# Create metadata dataframe
metadata <- merged_seurat@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Create sample column
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^NT1_"))] <- "NT1"
metadata$sample[which(str_detect(metadata$cells, "^NT2_"))] <- "NT2"
metadata$sample[which(str_detect(metadata$cells, "^PT1_"))] <- "PT1"
metadata$sample[which(str_detect(metadata$cells, "^PT2_"))] <- "PT2"
metadata$sample[which(str_detect(metadata$cells, "^PT3_"))] <- "PT3"
metadata$sample[which(str_detect(metadata$cells, "^RT1_"))] <- "RT1"
metadata$sample[which(str_detect(metadata$cells, "^RT2_"))] <- "RT2"
metadata$sample[which(str_detect(metadata$cells, "^RT3_"))] <- "RT3"

# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

# Create .RData object to load at any time
save(merged_seurat, file="data/merged_seurat_filtered.RData")


####### PLOTS ##########################
load("data/merged_seurat_filtered.RData")

#PT1 and RT3 are poor quality or too different to the other samples?

metadata<- merged_seurat@meta.data
##########################
#check data quality plotd
#https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html
# Visualize the number of cell counts per sample
metadata %>%
        ggplot(aes(x=sample, fill=sample)) +
        geom_bar() +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        theme(plot.title = element_text(hjust=0.5, face="bold")) +
        ggtitle("NCells")

#nUMI = nCount_RNA
#nGene = nFeature_RNA
# Visualize the number UMIs/transcripts per cell
metadata %>%
        ggplot(aes(color=sample, x=nCount_RNA, fill= sample)) +
        geom_density(alpha = 0.2) +
        scale_x_log10() +
        theme_classic() +
        ylab("Cell density") +
        geom_vline(xintercept = 500)


# Visualize the distribution of genes detected per cell via histogram
metadata %>%
        ggplot(aes(color=sample, x=nFeature_RNA, fill= sample)) +
        geom_density(alpha = 0.2) +
        theme_classic() +
        scale_x_log10() +
        geom_vline(xintercept = 300)

# Visualize the distribution of genes detected per cell via boxplot
metadata %>%
        ggplot(aes(x=sample, y=log10(nFeature_RNA), fill=sample)) +
        geom_boxplot() +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        theme(plot.title = element_text(hjust=0.5, face="bold")) +
        ggtitle("NCells vs NGenes")


# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>%
        ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) +
        geom_point() +
        scale_colour_gradient(low = "gray90", high = "black") +
        stat_smooth(method=lm) +
        scale_x_log10() +
        scale_y_log10() +
        theme_classic() +
        geom_vline(xintercept = 500) +
        geom_hline(yintercept = 250) +
        facet_wrap(~sample)

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>%
        ggplot(aes(color=sample, x=mitoRatio, fill=sample)) +
        geom_density(alpha = 0.2) +
        scale_x_log10() +
        theme_classic() +
        geom_vline(xintercept = 0.25)


# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
        ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
        geom_density(alpha = 0.2) +
        theme_classic() +
        geom_vline(xintercept = 0.8)





