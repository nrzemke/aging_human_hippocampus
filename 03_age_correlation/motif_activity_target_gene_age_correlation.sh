#!/bin/bash

#SBATCH -J motif_name #Optional, short for --job-name
#SBATCH -N 1 #Number of nodes
#SBATCH -n 1 #Total number of tasks increase this number to increase parallelization
#SBATCH -c 8 #Number of threads per process
#SBATCH --mem-per-cpu 8G
#SBATCH -t 1:00:00 #Short for --time walltime limit
#SBATCH -o /tscc/nfs/home/nzemke/cluster-logs/cellname_motif_name.out #standard output name
#SBATCH -e /tscc/nfs/home/nzemke/cluster-logs/cellname_motif_name.err #Optional, standard error name
#SBATCH -p condo #Partition name
#SBATCH -q condo #QOS name
#SBATCH -A csd788 #Allocation name

cd /tscc/projects/ps-renlab2/nzemke/multiome/hippocampus/40_donor_analysis/tf_fimo/

celltype=('Oligo' 'OPC' 'SUB' 'CA1' 'Macro' 'Endo' 'VLMC' 'CA2-CA3' 'SST' 'VIP' 'LAMP5' 'Chandelier' 'NR2F2' 'PVALB' 'T-Cell' 'Microglia' 'Astro' 'DG')

for celltype in ${celltype[@]}

do

mkdir motif_atac_pcc/${celltype}
mkdir motif_atac_pcc/${celltype}/rna_targets

motif_dir='/tscc/nfs/home/nzemke/renlab2/multiome/mop_m1_integration/human_atac/motif_epi_conservation/motif_files/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme/'

for x in ${motif_dir}* ; do
  motif=$(echo $x | sed "s=${motif_dir}==g");
  echo $motif


sed 1d fimo_motif_out/${motif}/fimo.tsv | cut -f3 | sed 's/:/\t/g' | sed 's/-/\t/g' | grep -Fw -f - ~/renlab2/multiome/hippocampus/40_donor_analysis/ATAC_peaks/Human.hippocampus.subclass.${celltype}_peaks.bed | cut -f1-4 > fimo_motif_out/${motif}/${celltype}_${motif}_peaks.bed

cut -f1-3 fimo_motif_out/${motif}/${celltype}_${motif}_peaks.bed | sed 's/\t/-/g' | grep -Fw -f - ~/renlab2/multiome/hippocampus/40_donor_analysis/correlation_with_age/ATAC_age_correlation/${celltype}_ATAC_pcc_donor_counts_filt_donors.tsv | sed 's/"//g' | cat header_atac - > motif_atac_pcc/${celltype}/${motif}_ATAC_pcc.tsv

cut -f1-3 fimo_motif_out/${motif}/${celltype}_${motif}_peaks.bed | grep -Fw -f - ~/renlab2/multiome/hippocampus/40_donor_analysis/abc/abc_predictions/${celltype}/EnhancerPredictions.txt | cut -f5 | sort | uniq | grep -Fw -f - ~/renlab2/multiome/hippocampus/40_donor_analysis/correlation_with_age/${celltype}_pcc_donor_20k_counts_2x_donors.tsv | cat header_rna - | sed 's/"//g' > motif_atac_pcc/${celltype}/rna_targets/${motif}_RNA_pcc.tsv

done

done
