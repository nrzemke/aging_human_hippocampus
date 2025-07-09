#!/bin/bash

motif_dir='/tscc/nfs/home/nzemke/renlab2/multiome/mop_m1_integration/human_atac/motif_epi_conservation/motif_files/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme/'

counter=0
for x in ${motif_dir}* ; do 
    if [ $counter -ge 1 ] && [ $counter -lt 792 ]; then
  motif=$(echo $x | sed "s=${motif_dir}==g");
echo $motif
   elif [ $counter -ge 792 ]; then
        break
    fi
    ((counter++))

fimo --o fimo_motif_out/${motif} ${motif_dir}${motif} union_peaks.fa

awk 'BEGIN{FS=OFS="\t"}($9<0.05){print $3,$4,$5,$6}' fimo_motif_out/${motif}/fimo.tsv | awk 'BEGIN {OFS=FS="\t"} {gsub("[:-]", "\t", $1); print}' | awk 'BEGIN{FS=OFS="\t"}{print $1,$2+$4,$2+$5,$6}' > fimo_motif_out/${motif}/motif_q.05.bed

done
