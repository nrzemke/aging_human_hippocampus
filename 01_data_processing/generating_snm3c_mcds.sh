# Generate mcds file for each sample in terminal
allcools generate-dataset \
    --allc_table \
    ${outfolder}/${samplename}/fastq_demultiplex_qc/AllcPaths.tsv \
    --output_path \
    ${outfolder}/${samplename}/allcools_analysis/${samplename}Cell.mcds \
    --chrom_size_path \
    ~/genome/hg38.nochrM.chrom.sizes \
    --obs_dim cell \
    --cpu 32 \
    --chunk_size 50 \
    --regions chrom100k 100000 \
    --regions chrom5k 5000 \
    --regions geneslop2k \
    ~/genome/hg38_gencode_v28.geneslop2k.sorted.bed.gz \
    --regions promoter \
    ~/genome/hg38_gencode_v28.promoter.sorted.bed.gz \
    --regions CGI \
    ~/genome/hg38_CGI.sorted.bed.gz \
    --quantifiers chrom100k count CGN,CHN,CAN \
    --quantifiers chrom5k count CGN,CHN,CAN \
    --quantifiers geneslop2k count CGN,CHN,CAN \
    --quantifiers promoter count CGN,CHN,CAN \
    --quantifiers CGI count CGN,CHN,CAN \
    --quantifiers chrom5k hypo-score CGN cutoff=0.9 \
    --quantifiers promoter hypo-score CGN cutoff=0.9 \
    --quantifiers CGI hypo-score CGN cutoff=0.9
