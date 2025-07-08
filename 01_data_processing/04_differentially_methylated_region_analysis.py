# # Differentially methylated region analysis
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

# allc paths for all samples
allc_paths = {}
for sample in samples:
    allc_path = '~/scmethylhic/human_hippocampus/snm3c/' + sample + '_deep/fastq_demultiplex/stats/AllcPaths.tsv'
    allc_paths[sample] = allc_path

metadata_path = 'metadata.csv.gz'
metadata = pd.read_csv(metadata_path, index_col=0)

clusters = ['Oligo', 'OPC', 'Astro', 'Micro', 'SUB', 'DG', 'CA', 'NR2F2-LAMP5', 'PVALB', 'SST', 'VIP', 'Endo_VLMC']
groups = ['2040', '4060', '6080', '80100']

for cluster in clusters:
    allcpath_c = metadata[metadata["cluster"] == cluster]
    
    for group in groups:
        allcpath_cg = allcpath_c[allcpath_c["group"] == group]
        
        # Define the output path
        output_path = p + 'mc_bulk/group/cluster/' + cluster + '_' + group + '.tsv'
        
        # Save 'AllcPath' values to file
        allcpath_cg[['AllcPath']].to_csv(output_path, sep="\t", index=False, header=False)

# in terminal
# 1. cluster
# 1) merge single allc files
conda activate allcools
clusters=('Oligo' 'OPC' 'Astro' 'Micro' 'SUB' 'DG' 'CA' 'NR2F2-LAMP5' 'PVALB' 'SST' 'VIP' 'Endo_VLMC')

for cluster in "${clusters[@]}"; do
    allcools merge-allc --allc_paths ${cluster}.tsv \
                        --output_path ${cluster}.tsv.gz \
                        --chrom_size_path /~/genome/hg38.nochrM.nochrUn.chrom.sizes
done

# 2) extract methylatioin context
for cluster in "${clusters[@]}"; do
     allcools extract-allc --allc_path ${cluster}.tsv.gz \
                           --output_prefix ${cluster} \
                           --mc_contexts CGN CHN \
                           --chrom_size_path ~/genome/hg38.nochrM.nochrUn.chrom.sizes
done

# 3) generate bigwig file from allc file
for cluster in "${clusters[@]}"; do
     allcools allc-to-bigwig --allc_path ${cluster}.tsv.gz \
                             --output_prefix ${cluster} \
                             --mc_contexts CGN CHN \
                             --chrom_size_path ~/genome/hg38.nochrM.nochrUn.chrom.sizes
done
# 2. cluster x age group
# 1) merge single allc files
conda activate allcools
clusters=('Oligo' 'OPC' 'Astro' 'Micro' 'SUB' 'DG' 'CA' 'NR2F2-LAMP5' 'PVALB' 'SST' 'VIP' 'Endo_VLMC')
groups=('2040' '4060' '6080' '80100')

for cluster in "${clusters[@]}"; do
  for group in "${groups[@]}"; do
    allcools merge-allc --allc_paths ${cluster}_${group}.tsv \
                        --output_path ${cluster}_${group}.tsv.gz \
                        --chrom_size_path /~/genome/hg38.nochrM.nochrUn.chrom.sizes
  done
done

# 2) extract methylatioin context
for cluster in "${clusters[@]}"; do
  for group in "${groups[@]}"; do
     allcools extract-allc --allc_path ${cluster}_${group}.tsv.gz \
                           --output_prefix ${cluster}_${group} \
                           --mc_contexts CGN CHN \
                           --chrom_size_path ~/genome/hg38.nochrM.nochrUn.chrom.sizes
  done
done

# 3) generate bigwig file from allc file
for cluster in "${clusters[@]}"; do
  for group in "${groups[@]}"; do
     allcools allc-to-bigwig --allc_path ${cluster}_${group}.tsv.gz \
                             --output_prefix ${cluster}_${group} \
                             --mc_contexts CGN CHN \
                             --chrom_size_path ~/genome/hg38.nochrM.nochrUn.chrom.sizes
  done
done

# in py
# 1. cluster
# Call differentially methylated sites & regions
import pandas as pd
import pathlib
from ALLCools.mcds import RegionDS
from ALLCools.dmr import call_dms, call_dmr, collapse_replicates

# Parameters
mc_bulk_dir = '~/scmethylhic/human_hippocampus/mc_bulk'
# make a dict, key is sample name, value is allc path
import pathlib
allc_table = {
    allc_path.name.split('.')[0]: str(allc_path) # change this based on your file names
    for allc_path in pathlib.Path(mc_bulk_dir).glob(
        '*.CGN-Both.allc.tsv.gz')
}
samples = list(allc_table.keys())
allc_paths = list(allc_table.values())

chrom_size_path = '~/genome/hg38.nochrM.nochrUn.chrom.sizes' 
output_dir = 'dmr' 

# Calculate CpG differentially methylated sites 
# Call differentially methylated sites
call_dms(
    output_dir=output_dir,
    allc_paths=allc_paths,
    samples=samples,
    chrom_size_path=chrom_size_path,
    cpu=45,
    max_row_count=50,
    n_permute=3000,
    min_pvalue=0.01)

# Call differentially methylated regions
call_dmr(output_dir=output_dir,
         p_value_cutoff=0.001,
         frac_delta_cutoff=0.3,
         max_dist=500,  
         residual_quantile=0.7,
         corr_cutoff=0.3,
         cpu=45
        )

# Open DMR matrix
dmr_dr = RegionDS.open(output_dir, region_dim='dmr') 
dmr_dr

# Annotate region DS
# DMR Overlapping Genome Features (BED)
# overlap the DMR regions with a set of BED files that containing different kinds of genome features (e.g. CGI, promoter)
genome_feature_dir = '~/genome/'
genome_feature_beds = {
    '.'.join(p.name.split('.')[:-3]): str(p) # Change this based on your file names
    for p in pathlib.Path(genome_feature_dir).glob('hg38_*.bed.gz')
}
beds = pd.Series(genome_feature_beds)
beds.to_csv('~/scmethylhic/human_hippocampus/mc_bulk/dmr/dmr_genome_featue_bed.csv', header=False)

dmr_dr.annotate_by_beds(slop=500,
                        bed_table='~/scmethylhic/human_hippocampus/mc_bulk/dmr/dmr_genome_featue_bed.csv',
                        dim='genome-features',
                        bed_sorted=False,
                        cpu=45)

# Filter DMRs
region_ds = RegionDS.open('dmr', select_dir=['dmr', 'dmr_genome-features'])

# Remove blacklist overlapping regions
is_blacklist= region_ds.get_feature('hg38_blacklist', 'genome-features').astype(bool)
# Filter by change DMR state
# change the blacklist-overlapping DMR state to 0 so they will not be included in the following selections
region_ds['dmr_state'].loc[{'dmr': is_blacklist.values}] = 0
new_region_ds = region_ds.sel({'dmr': ~is_blacklist.values})

# Get hypo- or hyper- state matrix
clusters = ['Oligo', 'OPC', 'Astro', 'Micro', 'SUB', 'DG', 'CA', 'NR2F2-LAMP5', 'PVALB', 'SST', 'VIP', 'Endo_VLMC']

# Initialize dictionaries to store hypo and hyper indices for each cluster
hypo_indices_dict = {}
hyper_indices_dict = {}

# Iterate over each cluster
for cluster in clusters:
    # Get hypo and hyper indices for the current cluster
    hypo_indices, hyper_indices = new_region_ds.get_hypo_hyper_index(cluster, use_collapsed=False)
    
    # Store hypo and hyper indices in the dictionaries
    hypo_indices_dict[cluster] = hypo_indices
    hyper_indices_dict[cluster] = hyper_indices

# 2. cluster x age
mc_bulk_group_dir = '~/scmethylhic/human_hippocampus/mc_bulk'

allc_table = {}
for allc_path in pathlib.Path(mc_bulk_group_dir).glob('*.CGN-Both.allc.tsv.gz'):
    sample_name = allc_path.name.split('.')[0]
    if sample_name in allc_table:
        allc_table[sample_name].append(str(allc_path))
    else:
        allc_table[sample_name] = [str(allc_path)]

samples = list(allc_table.keys())
allc_paths = [path for paths in allc_table.values() for path in paths]

chrom_size_path = '~/genome/hg38.nochrM.nochrUn.chrom.sizes' 
output_dir = '~/scmethylhic/human_hippocampusmc_bulk/dmr' 

# call DMSs
call_dms(
    output_dir=output_dir,
    allc_paths=allc_paths,
    samples=samples,
    chrom_size_path=chrom_size_path,
    cpu=45,
    max_row_count=50,
    n_permute=3000,
    min_pvalue=0.01)

# Call DMRs
call_dmr(output_dir=output_dir,
         p_value_cutoff=0.001,
         frac_delta_cutoff=0, 
         max_dist=500,  
         residual_quantile=0.7,
         corr_cutoff=0.3,
         cpu=45
        )

output_dir = 'dmr'
dmr_dr = RegionDS.open(output_dir, region_dim='dmr')
genome_feature_dir = '~/genome/'
genome_feature_beds = {
    '.'.join(p.name.split('.')[:-3]): str(p) # Change this based on your file names
    for p in pathlib.Path(genome_feature_dir).glob('hg38_*.bed.gz')
}
beds = pd.Series(genome_feature_beds)
beds.to_csv('~/scmethylhic/human_hippocampus/mc_bulk/dmr/dmr_genome_featue_bed.csv', header=False)

dmr_dr.annotate_by_beds(slop=500,
                        bed_table='~/scmethylhic/human_hippocampus/mc_bulk/dmr/dmr_genome_featue_bed.csv',
                        dim='genome-features',
                        bed_sorted=False,
                        cpu=45)

# Filter DMRs
region_ds = RegionDS.open('dmr', select_dir=['dmr', 'dmr_genome-features'])
is_blacklist= region_ds.get_feature('hg38_blacklist', 'genome-features').astype(bool)
region_ds['dmr_state'].loc[{'dmr': is_blacklist.values}] = 0
new_region_ds = region_ds.sel({'dmr': ~is_blacklist.values})

# Filter DMRs with less than 2 DMSs
result_df = pd.DataFrame(index=new_region_ds['dmr'].values)

# Initialize the 'hypo' and 'hyper' columns with empty strings
result_df['hypo'] = ''
result_df['hyper'] = ''

# Iterate over each cluster
for cluster in clusters:
    # Get hypo and hyper indices for the current cluster
    hypo_indices = hypo_indices_dict[cluster]
    hyper_indices = hyper_indices_dict[cluster]
    
    # Update the "hypo" column for DMRs identified as hypo-methylated in the current cluster
    result_df.loc[result_df.index.isin(hypo_indices), 'hypo'] += ',' + cluster
    
    # Update the "hyper" column for DMRs identified as hyper-methylated in the current cluster
    result_df.loc[result_df.index.isin(hyper_indices), 'hyper'] += ',' + cluster

# Remove leading comma from each entry in the 'hypo' and 'hyper' columns
result_df['hypo'] = result_df['hypo'].str.lstrip(',')
result_df['hyper'] = result_df['hyper'].str.lstrip(',')

data = {
    'dmr_chrom': new_region_ds['dmr_chrom'].values,
    'dmr_start': new_region_ds['dmr_start'].values,
    'dmr_end': new_region_ds['dmr_end'].values,
    'dmr_ndms': new_region_ds['dmr_ndms'].values,
    'dmr_length': new_region_ds['dmr_length'].values,
    'dmr_da_sum_mc': dmr_da_sum_mc,
    'dmr_da_sum_cov': dmr_da_sum_cov,
    'dmr':  new_region_ds['dmr'].values,
    'dmr_id': ['Union_DMR_' + str(i+1) for i in range(len(new_region_ds['dmr_chrom']))],
    'hypo': result_df['hypo'].values,
    'hyper': result_df['hyper'].values
}

df = pd.DataFrame(data)
condition = df['dmr_ndms'] >= 2
df_filtered = df[condition]
