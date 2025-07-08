# Integration wtih snRNA-seq data
import xarray as xr
import pandas as pd
import seaborn as sns
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import pybedtools
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

# Preprocessing mC dataset
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

# Create geneslop2k mcds 
# Create mcds paths for all samples
mcds_paths = [p + sample + '_deep/allcools_analysis/' + sample + 'Cell.mcds/' for sample in samples]

# Open MCDS files
mcds_list = []

for i in range(len(samples)):
    mcds = MCDS.open(
        mcds_paths[i],
        obs_dim='cell',
        var_dim='geneslop2k',
        use_obs=metadatas[i]['index']
    )
    mcds_list.append(mcds)

# before concat mcds, add sample information
mcds_list_con = []

for i in range(len(mcds_list)):
    sample = samples[i]
    mcds = mcds_list[i].copy()
    mcds['cell'] = sample + '-' + mcds.indexes['cell'] 
    mcds_list_con.append(mcds)

# Concatenate mcds
mcds = xr.concat(mcds_list_con, 'cell')
mcds.attrs = {"obs_dim":"cell", "var_dim":"geneslop2k"}

# Calculate gene mCG fraction
chrom_to_remove = ['chrM']

metadata_path = 'metadata.csv.gz'
gene_meta_path = '~/scmethylhic/human_hippocampus/GeneMetadata.csv.gz'

obs_dim = 'cell'
var_dim = 'geneslop2k'
min_cov = 5

# Load meta data
metadata = pd.read_csv(metadata_path, index_col=0)

# Load gene metadata
gene_meta = pd.read_csv(gene_meta_path, sep='\,', header=0, index_col=0)

# Filter genes by overlap and chromosomes
genes_to_skip = set()
# skip smaller genes mostly covered by a larger gene, e.g., a miRNA within a protein coding gene.
# F=0.9 means > 90% of gene_b is overlapped with gene_a, in this case, we only keep gene_a for DMG test
gene_bed = pybedtools.BedTool.from_dataframe(
    gene_meta.reset_index()[['chrom', 'start', 'end', 'gene_id']])
mapped_bam = gene_bed.map(b=gene_bed, c=4, o='distinct', F=0.9)
for _, (*_, gene_a, gene_b_str) in mapped_bam.to_dataframe().iterrows():
    for gene_b in gene_b_str.split(','):
        if gene_b != gene_a:
            genes_to_skip.add(gene_b)
# remove certain chromosomes
genes_to_skip |= set(gene_meta.index[gene_meta['chrom'].isin(chrom_to_remove)])
use_features = gene_meta.index[~gene_meta.index.isin(genes_to_skip)]
print(f'{use_features.size} features remained after filtering')

# Filter genes by cell mean coverage
mcds.add_feature_cov_mean()

feature_cov_mean = mcds.coords[f'{var_dim}_cov_mean'].to_pandas()
use_features &= feature_cov_mean[feature_cov_mean > min_cov].index

print(f'{use_features.size} features with feature coverage mean > {min_cov} remained')

mcds.filter_feature_by_cov_mean(min_cov=min_cov)

# Calculate and Save Gene mC Fractions
mcds.add_mc_frac(normalize_per_cell=True, clip_norm_value=10)
mcds.add_cell_metadata(metadata)

mcds = mcds[['geneslop2k_da_frac']]
mcds['geneslop2k_da_frac'] = mcds['geneslop2k_da_frac'].astype('float32')
mcds.write_dataset('geneslop2k_da_frac.mcds', var_dims=['geneslop2k'])

# Select Highly Variable Features (HVF)
mcds.calculate_hvf(mc_type='CGN', 
                   var_dim='geneslop2k', 
                   min_mean=0,
                   max_mean=5,
                   bin_min_features=5,
                   n_top_feature=25000,
                   mean_binsize=0.05,
                   cov_binsize=100)

# Get cell-by-feature mC fraction AnnData
adata = mcds.get_adata(mc_type='CGN', var_dim='geneslop2k', select_hvf=True) 

selected_var_names = adata.var_names[adata.var_names.isin(gene_meta.index)]
adata_sub = adata[:, selected_var_names]
gene_names = gene_meta.loc[adata_sub.var_names, 'gene_name'].tolist()
adata_sub.var_names = gene_names

# scale
log_scale(adata_sub)

# PCA
col_name = 'snm3c_celltype'

sc.tl.pca(adata_sub)
n_components = significant_pc_test(adata_sub)
fig, axes = plot_decomp_scatters(adata_sub,
                                 n_components=n_components,
                                 hue=col_name,
                                 hue_quantile=(0.25, 0.75),
                                 nrows=3,
                                 ncols=5)

# Calculate Nearest Neighbors
knn = -1 # -1 means auto determine

if knn == -1:
    knn = max(15, int(np.log2(adata_sub.shape[0])*2))
sc.pp.neighbors(adata_sub, n_neighbors=knn)

# Leiden Clustering
resolution = 0.2
sc.tl.leiden(adata_sub, resolution=resolution)

# UMAP
sc.tl.umap(adata_sub)

fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
_ = categorical_scatter(data=adata_sub,
                        ax=ax,
                        coord_base='umap',
                        hue='leiden',
                        text_anno='leiden',
                        show_legend=True)
adata_sub.write_h5ad('allcools_geneslop2k_CGN.h5ad')

# Perform cross-modality integration
# Load snRNAseq datasets and select shared HVFs
tenx_n_genes = 30000

tenx_adata = anndata.read_h5ad(f'~/scmethylhic/human_hippocampus/rna/human_aging_rna_allcools.h5ad')
tenx_hvf = tenx_adata.var['dispersions_norm'].sort_values(ascending=False).dropna()[:tenx_n_genes].index

snm3c_adata = anndata.read_h5ad('allcools_geneslop2k_CGN.h5ad')

hvfs = tenx_hvf.intersection(snm3c_adata.var_names)
with open('hvgs.txt', 'w') as f:
    f.write('\n'.join(hvfs))
print(hvfs.size, 'genes occured in all datasets.')

# Apply HVF selection and scale features
tenx_adata = tenx_adata[:, hvfs].copy()
sc.pp.scale(tenx_adata)

snm3c_adata = snm3c_adata[:, hvfs].copy()
sc.pp.scale(snm3c_adata)
snm3c_adata.X *= -1

# Concatenate all datasets
adata = tenx_adata.concatenate([
    snm3c_adata,
], index_unique='-')
adata.obs['batch'] = adata.obs['batch'].map({
    '0': '10X',
    '1': 'snm3C',
})

# PCA and Harmony Correction
sc.pp.pca(adata, n_comps=100)
sc.pl.pca_variance_ratio(adata)
n_components = significant_pc_test(adata, p_cutoff=0.2)

ho = run_harmony(adata.obsm['X_pca'],
                 meta_data=adata.obs,
                 vars_use='batch',
                 random_state=0,
                 nclust=100,
                 max_iter_harmony=20)

adata.obsm['X_pca'] = ho.Z_corr.T

_ = plot_decomp_scatters(adata,
                         n_components=n_components,
                         hue='batch',
                         palette='tab10')

# Co-Clustering
sc.pp.neighbors(adata, n_neighbors=25)
sc.tl.leiden(adata, resolution=1.5)
sc.tl.umap(adata)

adata.obs['umap_0'] = adata.obsm['X_umap'][:, 0]
adata.obs['umap_1'] = adata.obsm['X_umap'][:, 1]
adata.write_h5ad('allcools_snm3c_rna_integrated.h5ad')

# Label transferring from RNA to mC
import time
import anndata
import scanpy as sc
import pandas as pd
import seaborn as sns
import numpy as np
import scanpy.external as sce
import matplotlib as mpl
import matplotlib.pyplot as plt
import pynndescent

from glob import glob
from scipy.sparse import csr_matrix
from ALLCools.plot import *
from ALLCools.clustering import *
from ALLCools.integration.seurat_class import SeuratIntegration
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import normalize, OneHotEncoder
from matplotlib.colors import LogNorm
from matplotlib import cm as cm

# Load Datasets and select shared HVFs
adata = anndata.read_h5ad(f'~/scmethylhic/human_hippocampus/allcools_snm3c_rna_integrated.h5ad')

rna_cell = (adata.obs['batch']=='10X')
mc_cell = (adata.obs['batch']=='snm3C')

start_time = time.time()
index = pynndescent.NNDescent(adata.obsm['X_pca'][rna_cell], metric='euclidean', n_neighbors=50, random_state=0, n_jobs=-1)
G, D = index.query(adata.obsm['X_pca'][mc_cell], k=15)
chunk_size = 50000
sd = 1
start_time = time.time()
cellfilter = D[:, -1] == 0
D = 1 - D / D[:, -1][:, None]
D[cellfilter] = 1
D = 1 - np.exp(-D * (sd**2) / 4)
D = D / (np.sum(D, axis=1) + 1e-6)[:, None]

rna_index = rna_cell.index[rna_cell]
mc_index = mc_cell.index[mc_cell]

enc = OneHotEncoder()
labelref = enc.fit_transform(adata.obs.loc[rna_index, ['Final1']].astype(str)).toarray()
cluster = pd.DataFrame(index=mc_index, columns=['rnatype', 'score'], dtype=str)

for chunk_start in range(0, len(mc_index), chunk_size):
    result = (
        D[chunk_start : (chunk_start + chunk_size), :, None]
        * labelref[G[chunk_start : (chunk_start + chunk_size)].flatten(), :]
        .reshape((-1, 15, enc.categories_[0].shape[0]))
    ).sum(axis=1)

    result = pd.DataFrame(
        result,
        columns=enc.categories_[0],
        index=mc_index[chunk_start : (chunk_start + chunk_size)],
    )

    result = result.loc[:, result.columns != "nan"]
    cluster.loc[result.index, "rnatype"] = result.idxmax(axis=1).values
    cluster.loc[result.index, "score"] = result.max(axis=1).values
    print(chunk_start)

print(time.time() - start_time)
cluster.to_hdf("mc_rnacluster.hdf", key="data")

cluster['dipctype'] = adata.obs.loc[mc_index, 'snm3c_celltype']
adata.obs.loc[mc_index, 'cell-type'] = cluster['rnatype']

adata_snm3c = adata[adata.obs['batch'] == 'snm3C'].copy()
