#########################################################################################
### Seurat workflow for all cells.		
### The FASTQ files were mapped to the GRCh38 human reference genome and mm10 mouse reference genome to distinguish human and mouse cell using CellRanger count.
### Control and Gem group were aggregated using the cellranger aggr pipeline.
#########################################################################################

library(Seurat)
# Load aggregated data
gem.data <- Read10X(data.dir = "filtered_feature_bc_matrix/")
gem <- CreateSeuratObject(counts = gem.data, project = "PDX_10X", min.cells = 3, min.features = 300)

# Add group information
gem@meta.data$group=NA
gem@meta.data$group[grep("pdx_gem_",rownames(gem))]="gem"
gem@meta.data$group[grep("pdx_con_",rownames(gem))]="con"
sce <- gem

# QC and selecting cells for further analysis
sce <- PercentageFeatureSet(sce, pattern = "^Mt-", col.name = "mitoRatio")
VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), ncol = 3)
sce <- subset(sce, subset = nFeature_RNA > 300 & nCount_RNA >1000 & mitoRatio<0.2 &  log10GenesPerUMI>0.8  )

# Apply SCTransform normalization
sce <- SCTransform(sce, vars.to.regress = "mitoRatio", verbose = FALSE)

# Perform dimensionality reduction by PCA and UMAP embedding
sce <- RunPCA(sce)
sce <- RunUMAP(sce, dims = 1:40)
sce <- FindNeighbors(sce, dims = 1:40)
sce <- FindClusters(sce)

# Visualization
DimPlot(sce, reduction = "tsne", label = TRUE)
DimPlot(sce, reduction = "tsne", label = TRUE, split.by = "sample")

# Save data
save(sce, file = "10x_PDX_sce.tsne.rds")