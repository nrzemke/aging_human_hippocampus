#!/bin/bash

## Convert HiC from mcool to bedpe format for ABC compatibility 
## source for python script: https://github.com/zsq-berry/3D-genome-tools/blob/main/cool2hic.py


#conda activate ~/renlab2/conda_envs/cooler

for x in *.cool; do
  type=$(echo $x | sed 's=.Q.cool==g');
  echo ${type}

mkdir ${type}

python3 cool2hic.py -i ${type}.Q.cool -r 10000


{

chrom=('1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22')

for chrom in ${chrom[@]}

do

mkdir ${type}/chr${chrom}

grep -w "chr${chrom}" matrix.txt | awk 'BEGIN{FS=OFS="\t"}{print $2,$3,$3+10000,$6,$7,$7+10000,".",$9}' > ${type}/chr${chrom}/chr${chrom}.bedpe

gzip ${type}/chr${chrom}/chr${chrom}.bedpe


done

}

rm matrix.txt

done

=====================================

## Call neighborhoods

#!/bin/bash

#module load samtools
#module load bedtools

for x in ../ATAC_bam/*_sorted.bam ; do 
  type=$(echo $x | sed 's=../ATAC_bam/==g' | sed 's=_sorted.bam==g');
  echo ${type}


python ~/renlab2/software/ABC-Enhancer-Gene-Prediction/src/run.neighborhoods.py \
--candidate_enhancer_regions ~/renlab2/multiome/hippocampus/40_donor_analysis/ATAC_peaks/Human.hippocampus.subclass.${type}_peaks.bed \
--genes hg38_gene_bounds_protein_coding.txt \
--ATAC ../ATAC_bam/${type}_sorted.bam \
--expression_table ~/renlab2/multiome/hippocampus/40_donor_analysis/rna_cpm/${type}_cpm.tsv \
--chrom_sizes hg38.chrom.sizes \
--ubiquitously_expressed_genes ubiqiutously_expressed_empty.txt \
--cellType ${type} \
--outdir neighborhoods/${type}

done

=====================================

## Make ABC predictions

#!/bin/bash

for x in mcool_to_bedpe/* ; do 
  type=$(echo $x);
  echo ${type}


mkdir abc_predictions/${type}

python ~/renlab2/software/ABC-Enhancer-Gene-Prediction/src/predict.py \
--enhancers neighborhoods/${type}/EnhancerList.txt \
--genes neighborhoods/${type}/GeneListfilt.txt \
--HiCdir mcool_to_bedpe/${type} \
--hic_type bedpe \
--chrom_sizes hg38.chrom.sizes.autosomes \
--hic_resolution 10000 \
--threshold .02 \
--cellType ${type} \
--outdir abc_predictions/${type} \
--make_all_putative

done

