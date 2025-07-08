# Load libraries
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)

# Optional memory limits for large datasets
options(future.globals.maxSize = 80000 * 1024^2)
options(Seurat.object.assay.version = "v5")

# ----------------------------
# PART 1: Data Preprocessing
# ----------------------------
# Split object by sample/donor
sample_list <- SplitObject(object, split.by = "orig.ident")

# Run SCTransform on each sample separately with scaling and centering disabled
sample_list <- lapply(sample_list, function(x) {
  SCTransform(x, do.scale = FALSE, do.center = FALSE, verbose = FALSE)
})

# Find 3000 variable genes per sample
sample_list <- lapply(sample_list, function(x) {
  FindVariableFeatures(x, nfeatures = 3000)
})

# Merge all samples into one object
object <- merge(sample_list[[1]], y = sample_list[-1])

# Get all variable genes from all samples
all_variable_genes <- unique(unlist(lapply(sample_list, function(x) {
  VariableFeatures(x)
})))

# Filter to keep only protein-coding genes
# You'll need to provide a list of protein-coding genes or use a reference
# Example using a hypothetical protein_coding_genes vector:
# protein_coding_genes <- read.csv("protein_coding_genes.csv")$gene_name
# protein_coding_variable_genes <- intersect(all_variable_genes, protein_coding_genes)

# For now, assuming all_variable_genes are protein-coding or you filter them:
protein_coding_variable_genes <- all_variable_genes  # Replace with actual filtering

# Select top 2000 protein-coding variable genes based on variance
# Get variance info for ranking
variance_info <- HVFInfo(object)
protein_coding_variance <- variance_info[protein_coding_variable_genes, ]
protein_coding_variance <- protein_coding_variance[order(protein_coding_variance$variance.standardized, decreasing = TRUE), ]
top_2000_genes <- rownames(protein_coding_variance)[1:min(2000, nrow(protein_coding_variance))]

# Set the top 2000 protein-coding genes as variable features
VariableFeatures(object) <- top_2000_genes

# Scale data and run PCA using the top 2000 protein-coding variable genes
object <- ScaleData(object)
object <- RunPCA(object)

# ----------------------------
# PART 2: rPCA Integration for Batch Correction
# ----------------------------
# Batch correction across donors using reciprocal PCA (rPCA) on SCTransformed PCs
object <- IntegrateLayers(
  object = object, 
  method = RPCAIntegration,
  orig.reduction = "pca", 
  new.reduction = "integrated.rpca",
  verbose = FALSE
)

# Build k-nearest neighbors graph using PCs 1:30 from rPCA-corrected space
object <- FindNeighbors(object, reduction = "integrated.rpca", dims = 1:30)

# Identify clusters using Leiden clustering
object <- FindClusters(object, resolution = 0.3, cluster.name = "rpca_clusters", algorithm = 4)  # algorithm = 4 for Leiden

# Run UMAP on rPCA-corrected space for visualization
object <- RunUMAP(object, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")

# ----------------------------
# Visualization
# ----------------------------
p1 <- DimPlot(
  object,
  reduction = "umap.rpca",
  group.by = c("orig.ident", "celltypes", "rpca_clusters"),
  combine = FALSE, 
  label.size = 2
)

# Display plots
p1

# ----------------------------
# PART 3: Cell Type Subsetting and Harmony Integration
# ----------------------------
# Example: Subset specific cell types for sub-clustering
# Replace "your_celltype" with actual cell type names
celltype_subset <- subset(object, subset = celltypes == "your_celltype")

# Re-normalize using SCTransform after subsetting
celltype_subset <- SCTransform(celltype_subset, verbose = FALSE)

# Standard preprocessing for the subset
celltype_subset <- FindVariableFeatures(celltype_subset)
celltype_subset <- ScaleData(celltype_subset)
celltype_subset <- RunPCA(celltype_subset)

# Apply Harmony for batch correction within cell type subset
celltype_subset <- IntegrateLayers(
  object = celltype_subset, 
  method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "harmony",
  verbose = FALSE
)

# Build k-nearest neighbors graph using Harmony-corrected PCs
celltype_subset <- FindNeighbors(celltype_subset, reduction = "harmony", dims = 1:30)

# Sub-clustering within cell type using Leiden clustering
celltype_subset <- FindClusters(celltype_subset, resolution = 0.5, cluster.name = "harmony_subclusters", algorithm = 4)

# Run UMAP on Harmony-corrected space for visualization
celltype_subset <- RunUMAP(celltype_subset, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

# Visualization of sub-clusters
p2 <- DimPlot(
  celltype_subset,
  reduction = "umap.harmony",
  group.by = c("orig.ident", "harmony_subclusters"),
  combine = FALSE, 
  label.size = 2
)

# Display plots
p2

# ----------------------------
# Alternative: Loop through multiple cell types
# ----------------------------
# If you want to process multiple cell types:
# unique_celltypes <- unique(object$celltypes)
# celltype_objects <- list()
# 
# for(celltype in unique_celltypes) {
#   # Subset cell type
#   temp_obj <- subset(object, subset = celltypes == celltype)
#   
#   # Re-normalize using SCTransform
#   temp_obj <- SCTransform(temp_obj, verbose = FALSE)
#   temp_obj <- FindVariableFeatures(temp_obj)
#   temp_obj <- ScaleData(temp_obj)
#   temp_obj <- RunPCA(temp_obj)
#   
#   # Apply Harmony for batch correction
#   temp_obj <- IntegrateLayers(
#     object = temp_obj, 
#     method = HarmonyIntegration,
#     orig.reduction = "pca", 
#     new.reduction = "harmony",
#     verbose = FALSE
#   )
#   
#   # Clustering and UMAP
#   temp_obj <- FindNeighbors(temp_obj, reduction = "harmony", dims = 1:30)
#   temp_obj <- FindClusters(temp_obj, resolution = 0.5, cluster.name = "harmony_subclusters", algorithm = 4)
#   temp_obj <- RunUMAP(temp_obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
#   
#   # Store in list
#   celltype_objects[[celltype]] <- temp_obj
# }

# ----------------------------
# Optional: Save final objects
# ----------------------------
# save(object, file = "main_seurat_object_with_rpca.RData")
# save(celltype_subset, file = "celltype_subset_with_harmony.RData")
