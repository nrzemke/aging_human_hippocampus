#!/bin/bash

# Author: Sainath Mamde
# Project: Cell-type and Age-Resolved ATAC-seq Peak Calling and Merging
# Genome: hg38
# Tools: MACS3, bedtools, awk, grep, sort

# Set working directory
WORKDIR="path/work/directory"
BLACKLIST="path/work/directory/hg38.blacklist.bed.gz"
CHROMSIZES="path/work/directory/hg38.chrom.sizes"

mkdir -p "$WORKDIR"
cd "$WORKDIR" || exit

# ----------------------------------------------
# Step 1: Peak Calling with MACS3 for BED Inputs
# ----------------------------------------------
echo "Step 1: Calling peaks with MACS3..."
for bed_file in *.bed; do
    celltype="${bed_file%.bed}"
    macs3 callpeak \
        --treatment "$bed_file" \
        --ext 150 --shift -75 --nomodel \
        -g hs \
        --name "$celltype" \
        -q 0.1 \
        --call-summits \
        -f BED \
        --outdir "$WORKDIR"
done

# ----------------------------------------------------
# Step 2: Filter Summits (Remove _alt and non-standard)
# ----------------------------------------------------
echo "Step 2: Filtering summits (remove _alt)..."
for file in *_summits.bed; do
    grep 'chr' "$file" | grep -v '_alt' > "${file%_summits.bed}_filtered_summits.bed"
done

# ----------------------------------------------------
# Step 3: Generate Celltype-Path Map for Merging Input
# ----------------------------------------------------
echo "Step 3: Creating input path map for merging..."
> celltype_paths.txt
for file in *_filtered_summits.bed; do
    celltype=$(basename "$file" _filtered_summits.bed)
    echo -e "$celltype\t$(realpath "$file")" >> celltype_paths.txt
done

# ----------------------------------------------------
# Step 4: Merge Filtered Peaks Using Iterative Overlap
# Found here: https://github.com/yal054/snATACutils/blob/2b62147e1378a6a1935ebc7fb1129320f4eb8754/bin/iterative_overlap_peak_merging.R#L2
# ----------------------------------------------------
echo "Step 4: Running iterative peak merging..."
Rscript iterative_overlap_peak_merging_1.R \
    -i "$WORKDIR/celltype_paths.txt" \
    -g hg38 \
    --blacklist "$BLACKLIST" \
    --chromSize "$CHROMSIZES" \
    -d "$WORKDIR" \
    -o all_liver_peaks
