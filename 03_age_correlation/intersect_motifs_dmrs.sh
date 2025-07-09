#!/bin/bash

celltype=('Oligo' 'OPC' 'SUB' 'CA' 'Endo_VLMC' 'SST' 'VIP' 'NR2F2-LAMP5' 'PVALB' 'Microglia' 'Astro' 'DG')

for celltype in ${celltype[@]}

do

mkdir motif_dmrs/${celltype}

motif_dir='/tscc/nfs/home/nzemke/renlab2/multiome/mop_m1_integration/human_atac/motif_epi_conservation/motif_files/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme/'

counter=0
for x in ${motif_dir}* ; do 
  motif=$(echo $x | sed "s=${motif_dir}==g");


sed 1d ../../correlation_with_age/dmr/${celltype}_dmr_mc_ratio_pearson.tsv | cut -f1,2 | sed 's/:/\t/' | sed 's/-/\t/' > ${celltype}_dmr_mc_ratio_pearson.bed

awk 'BEGIN{FS=OFS="\t"}{print $3,$4,$5,$6}' ../fimo_motif_out/${motif}/fimo.tsv | awk 'BEGIN {OFS=FS="\t"} {gsub("[:-]", "\t", $1); print}' | awk 'BEGIN{FS=OFS="\t"}{print $1,$2+$4,$2+$5}' | grep ^chr | bedtools intersect -u -a ${celltype}_dmr_mc_ratio_pearson.bed -b stdin > motif_dmrs/${celltype}/${motif}_pcc.tsv

done

done
