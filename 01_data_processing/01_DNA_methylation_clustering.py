# Create concatenated mcds and cluster 
import xarray as xr
import pandas as pd
import seaborn as sns
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import anndata
import bbknn

from ALLCools.mcds import MCDS
from ALLCools.clustering import tsne, significant_pc_test, log_scale
from ALLCools.plot import *
from ALLCools.clustering import \
    tsne, \
    significant_pc_test, \
    filter_regions, \
    remove_black_list_region, \
    lsi, \
    log_scale, \
    binarize_matrix
from annoy import AnnoyIndex

samples = [
"hc1153",
"hc11",
"hc1271",
"hc12",
"hc19",
"hc26",
"hc29",
"hc35",
"hc40",
"hc73",
"hc76_1",
"hc76_2",
"hc78",
"hc81",
"hc8",
"hc98",
"hc9",
"hs46426",
"k1203",
"k1216_1",
"k1216_2",
"k1265",
"m1134",
"m1745",
"m4781",
"m5021",
"m5087",
"M5276",
"m5551",
"m5579",
"m5610",
"m5614",
"m6052",
"m935",
"m937_1",
"m937_2",
"ms212191",
"ms69984_1",
"ms69984_2",
"ms73787",
"pb13344",
"pb13394",
"PB13414",
"s6021"
]

# Create metadata paths for all samples
# base path
p = '~/scmethylhic/human_hippocampus/snm3c/'

metadata_paths = []
for sample in samples:
    metadata_path = p + sample + '_deep/fastq_demultiplex_qc/metadata.passqc.doublet.csv.gz'
    metadata_paths.append(metadata_path)

metadatas = []
for path in metadata_paths:
    metadata = pd.read_csv(path, index_col=0)
    metadatas.append(metadata)

# Concatenate metadata
metadata_con = pd.concat(metadatas, axis=0)

# Create mcds paths for all samples
p = '~/scmethylhic/human_hippocampus/snm3c/'

mcds_paths = [p + sample + '_deep/allcools_analysis/' + sample + 'Cell.mcds/' for sample in samples]

mcds_list = []
for i in range(len(samples)):
    mcds = MCDS.open(
        mcds_paths[i],
        obs_dim='cell',
        var_dim='chrom5k',
        use_obs=metadatas[i]['index']
    )
    mcds_list.append(mcds)

# Before concat mcds, add sample information
mcds_list_con = []
for i in range(len(mcds_list)):
    sample = samples[i]
    mcds = mcds_list[i].copy()
    mcds['cell'] = sample + '-' + mcds.indexes['cell'] 
    mcds_list_con.append(mcds)

# Concatenate mcds 
mcds = xr.concat(mcds_list_con, 'cell')
mcds.attrs = {"obs_dim":"cell", "var_dim":"chrom5k"}

mcds.add_cell_metadata(metadata_reindex)
mcds

# Load hypo-methylation score matrix from MCDS
mcad = mcds.get_score_adata(mc_type='CGN', quant_type='hypo-score')

# Clustering based on mCG without batch effect correaction
# Binarize
binarize_matrix(mcad, cutoff=0.95)

# Filter Features
filter_regions(mcad, hypo_percent=0.1) # min % of cells that are non-zero in this region 

# TF-IDF Transform and Dimension Reduction
lsi(mcad, algorithm='arpack', obsm='X_pca', scale_factor=100000, n_components=50, random_state=0, fit_size=None)

# Choose significant components
pc_cutoff=0.05
n_components = significant_pc_test(mcad, p_cutoff=pc_cutoff, update=True)

hue = 'mCGFrac'
if hue in metadata.columns:
    mcad.obs[hue] = mcad.obs[hue].reindex(mcad.obs_names)
    fig, axes = plot_decomp_scatters(mcad,
                                     n_components=n_components,
                                     hue=hue,
                                     hue_quantile=(0.25, 0.75),
                                     nrows=5,
                                     ncols=5)

# Calculate Nearest Neighbors
knn = -1 # -1 means auto determine

if knn == -1:
    knn = max(15, int(np.log2(mcad.shape[0])*2))
    
sc.pp.neighbors(mcad, n_neighbors=knn)

# Leiden Clustering
resolution = 1
sc.tl.leiden(mcad, resolution=resolution)

# UMAP
sc.tl.umap(mcad)

# Clustering based on mCG with BBKNN batch effect correaction
%%time
sc.external.pp.bbknn(mcad, batch_key='sample', use_rep='X_pca', 
                     neighbors_within_batch=3, trim=100, # higher number, stronger batch correction
                     n_pcs=45) # n_pcs default=50 and this is too many for sc.tl.umap function here 

sc.tl.umap(mcad, min_dist=0.5, spread=1.5) 

# Subcluster neurons with mCH
subset_group = mcad.obs['leiden'].isin(['9','12','13','14','16'])
adata = mcad[subset_group, :]

# Binarize
binarize_matrix(adata, cutoff=0.95)

# Filter Features
filter_regions(adata)

# TF-IDF Transform and Dimension Reduction
lsi(adata, algorithm='arpack', obsm='X_pca')

# Choose significant components
pc_cutoff=0.05
significant_pc_test(adata, p_cutoff=pc_cutoff, update=True)

# Calculate Nearest Neighbors
knn = -1 # -1 means auto determine

if knn == -1:
    knn = max(15, int(np.log2(adata.shape[0])*2))
sc.pp.neighbors(adata, n_neighbors=knn)

# Leiden Clustering
resolution = 0.5
sc.tl.leiden(adata, resolution=resolution)

# UMAP
sc.tl.umap(adata)

%%time
sc.external.pp.bbknn(adata, batch_key='sample', use_rep='X_pca', n_pcs=25) 

sc.tl.umap(adata)

# Leiden Clustering
resolution = 1.1
sc.tl.leiden(adata, resolution=resolution)
