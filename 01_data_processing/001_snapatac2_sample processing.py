import snapatac2 as snap
import scanpy as sc

import os

# Define input fragment files and sample names
files = [
    ("Sample1", "path/to/sample1_fragments.tsv.gz"),
    ("Sample2", "path/to/sample2_fragments.tsv.gz"),
    # Add more samples as needed
]

# Step 1: Import fragments from multiple samples
adata = snap.pp.import_fragments(
    [fl for _, fl in files],
    file=[name + '.h5ad' for name, _ in files],
    chrom_sizes=snap.genome.hg38,
    min_num_fragments=500  # adjust if needed
)

# Step 2: Cell filtering
snap.pp.filter_cells(
    adata,
    min_counts=500,
    min_tsse=5,
    max_counts=1000000
)

# Step 3: Create tile matrix and feature selection
snap.pp.add_tile_matrix(adata)
snap.pp.select_features(adata, n_features=250000)

# Step 4: Doublet detection and removal
snap.pp.scrublet(adata)
snap.pp.filter_doublets(adata)

# Step 5: Dimensionality reduction and embedding
snap.tl.spectral(adata)
snap.tl.umap(adata)

# Optional: Clustering
sc.pp.neighbors(adata, use_rep="X_spectral")
sc.tl.leiden(adata, key_added="leiden")


adata.write("combined_cleaned_spectral.h5ad")
