# Load required libraries
library(Seurat)
library(DoubletFinder)
library(ggplot2)
library(tidyverse)

# Read 10X data (replace with your own path)
counts <- Read10X("path/to/raw_feature_bc_matrix")

# Create Seurat object
hippo <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200, project = "SampleProject")

# Calculate QC metrics
hippo[["percent.mt"]] <- PercentageFeatureSet(hippo, pattern = "^MT-")
hippo <- PercentageFeatureSet(hippo, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", col.name = "percent.ribo")

# Plot QC metrics
pdf("QC_metrics.pdf")
VlnPlot(hippo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 3)
dev.off()

# QC filtering
hippo <- subset(hippo, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 15)
hippo <- hippo[!grepl("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", rownames(hippo)), ]
hippo <- hippo[!grepl("^MT-", rownames(hippo)), ]

# Normalize and identify variable features
hippo <- NormalizeData(hippo)
hippo <- FindVariableFeatures(hippo, selection.method = "vst", nfeatures = 5000)

# Scaling and dimensionality reduction
hippo <- ScaleData(hippo)
hippo <- RunPCA(hippo, features = VariableFeatures(hippo))
hippo <- FindNeighbors(hippo, dims = 1:20)
hippo <- FindClusters(hippo, resolution = 0.5)
hippo <- RunUMAP(hippo, dims = 1:20)

# DoubletFinder: parameter sweep
sweep.res.list <- paramSweep_v3(hippo, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Plot pK selection
pdf("DoubletFinder_paramSweep.pdf")
ggplot(bcmvn, aes(pK, BCmetric, group = 1)) + geom_point() + geom_line()
dev.off()

# Optimal pK selection
pK <- bcmvn %>%
  filter(BCmetric == max(BCmetric)) %>%
  pull(pK) %>%
  as.numeric()

# Estimate expected doublet rate
annotations <- hippo@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.1 * nrow(hippo@meta.data))  # Adjust this expected rate per dataset
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

# Run DoubletFinder
hippo <- doubletFinder_v3(hippo, PCs = 1:20, pN = 0.25, pK = pK,
                         nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)

# Plot doublets
pdf("Doublets_umap.pdf")
DimPlot(hippo, reduction = "umap", group.by = colnames(hippo@meta.data)[ncol(hippo@meta.data)])
dev.off()

# Save processed object
save(hippo, file = "hippo_processed.RData")
