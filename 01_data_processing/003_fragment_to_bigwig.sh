#!/bin/bash

for x in /tscc/projects/ps-renlab2/smamde/Doublet_Corrected_TSCC/frag_peaks_cell_type_updated/*.bed ; do 
  type=$(echo $x | sed 's=/tscc/projects/ps-renlab2/smamde/Doublet_Corrected_TSCC/frag_peaks_cell_type_updated/==g' | sed 's=.bed==g');
  echo ${type}

/tscc/projects/ps-renlab2/nzemke/software/bedtools2/bin/bedToBam -i /tscc/projects/ps-renlab2/smamde/Doublet_Corrected_TSCC/frag_peaks_cell_type_updated/${type}.bed -g ~/renlab2/genome_files/hg38.chrom.sizes > atac_bam/${type}.bam

/tscc/projects/ps-renlab2/nzemke/software/samtools-1.19.2/samtools sort atac_bam/${type}.bam > atac_bam/${type}_sorted.bam

rm atac_bam/${type}.bam

/tscc/projects/ps-renlab2/nzemke/software/samtools-1.19.2/samtools index atac_bam/${type}_sorted.bam

bamCoverage --normalizeUsing RPKM -b atac_bam/${type}_sorted.bam -o bigwig/${type}.bw

done
