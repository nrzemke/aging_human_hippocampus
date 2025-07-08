# Differentially methylated gene analysis
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

# Identify differentially methylated genes - mCG
def get_gene_values_by_name(gene_name):
    data = gene_frac_da.sel(geneslop2k=gene_name_to_gene_id[gene_name]).to_pandas()
    data.name = gene_name
    return data

mcds_paths = 'geneslop2k_frac.mcds'
metadata_path = 'metadata.csv.gz'
gene_meta_path = '~/scmethylhic/human_hippocampus/GeneMetadata.csv.gz'

cluster_col = 'cluster'

obs_dim = 'cell'
var_dim = 'geneslop2k'
mc_type = 'CGN'  

top_n = 1000
auroc_cutoff = 0.8  # AUROC cutoff to report significant DMG
adj_p_cutoff = 0.05  # adjusted P value cutoff to report significant DMG
fc_cutoff = 0.8  # mC fraction fold change cutoff to report significant DMG
max_cluster_cells = 5000  # The maximum number of cells from a group, downsample large group to this number
max_other_fold = 2  # The fold of other cell numbers comparing
cpu = 30

# Load
metadata = pd.read_csv(metadata_path, index_col=0)

# Calculate DMG
dmg_table = one_vs_rest_dmg(metadata,
                            group=cluster_col,
                            mcds_paths=mcds_paths,
                            obs_dim=obs_dim,
                            var_dim=var_dim,
                            mc_type=mc_type,
                            top_n=top_n,
                            adj_p_cutoff=adj_p_cutoff,
                            fc_cutoff=fc_cutoff,
                            auroc_cutoff=auroc_cutoff,
                            max_cluster_cells=max_cluster_cells,
                            max_other_fold=max_other_fold,
                            cpu=cpu,
                            verbose=True)

dmg_table.to_hdf(f'{cluster_col}.OneVsRestDMG_CGN.hdf', key='data')
