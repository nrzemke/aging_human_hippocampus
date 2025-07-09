library(ggplot2)
library(pheatmap)
library(plyr)
library(viridis)
library(tidyr)
library(dplyr)
library(ggrepel)

cells <- c('Oligo', 'OPC', 'SUB', 'CA1', 'Macro', 'Endo', 'VLMC', 'CA2-CA3', 'SST', 'VIP', 'LAMP5', 'Chandelier', 'NR2F2', 'PVALB', 'T-Cell', 'Microglia', 'Astro', 'DG')

#Background sampling takes equal ratio of promoter proximal and promoter distal

prox <- read.table('~/renlab2/multiome/hippocampus/40_donor_analysis/ATAC_peaks/Human.hippocampus.all_peaks_all_cell_types_tss_prox_1kb.bed')
distal <- read.table('~/renlab2/multiome/hippocampus/40_donor_analysis/ATAC_peaks/Human.hippocampus.all_peaks_all_cell_types_tss_distal_1kb.bed')
prox$peak <- paste(prox[,1], prox[,2], prox[,3], sep = "-")
distal$peak <- paste(distal[,1], distal[,2], distal[,3], sep = "-")

prox[,1:4] <- NULL
distal[,1:4] <- NULL
prox$V5 <- 'Proximal'
distal$V5 <- 'Distal'
colnames(prox) <- c('group', 'peak')
colnames(distal) <- c('group', 'peak')
anno <- rbind(prox, distal)

for (celltype in cells) {

all <- read.table(paste0('~/renlab2/multiome/hippocampus/40_donor_analysis/correlation_with_age/ATAC_age_correlation/', celltype, '_ATAC_pcc_donor_counts_filt_donors.tsv'), header=T)
all$scaled <- scale(all$cor)
null <- all
null$peak <- row.names(null)
null <- join(null, anno, by = 'peak')
    
all_motifs <- read.table(paste0('motif_atac_pcc/', celltype, '/', celltype, '_motifs_100_ccres.tsv'),header=F)
motifs <- as.list(all_motifs$V1)

vol_val <- data.frame(fold = numeric(),
                      p_value = numeric(),
                      motif = character(),
                      stringsAsFactors = FALSE)

for (var in motifs) {
    
pcc <- read.table(paste0('motif_atac_pcc/', celltype, '/', var, '_ATAC_pcc.tsv'), header=T)

pcc <- join(pcc, null, by = 'peak')

ratio_proximal <- sum(pcc$group == "Proximal") / nrow(pcc)

# Calculate the number of rows to sample
num_proximal <- round(ratio_proximal * nrow(prox))
num_distal <- round((1 - ratio_proximal) * nrow(prox))

# Randomly sample rows from 'Proximal' group
proximal_sample <- null %>% 
  filter(group == "Proximal") %>% 
  sample_n(min(num_proximal, sum(null$group == "Proximal")), replace = FALSE)

# Randomly sample rows from 'Distal' group
distal_sample <- null %>% 
  filter(group == "Distal") %>% 
  sample_n(min(num_distal, sum(null$group == "Distal")), replace = FALSE)

# Combine the sampled data into a new dataframe
null_sampled <- bind_rows(proximal_sample, distal_sample)

    
vol_val <- rbind(vol_val, data.frame(scaled = mean(pcc$scaled, na.rm = TRUE),
                                      p_value = wilcox.test(null_sampled$scaled, pcc$scaled, na.rm = TRUE)$p.value,
                                      motif = var))
    
    }

vol_val$FDR <- p.adjust(vol_val$p_value, method = 'fdr')

options(repr.plot.width = 8, repr.plot.height = 8, repr.plot.res = 200)

vol_val$color <- ifelse(vol_val$scaled > 0, "#FC8961",
                         ifelse(vol_val$scaled < 0, "#7F2582", "grey"))

pos <- subset(vol_val, subset = scaled > 0)
top_positive <- head(pos[order(pos$FDR), ],10)
neg <- subset(vol_val, subset = scaled < 0)
top_negative <- head(neg[order(neg$FDR), ],10)
top_20 <- rbind(top_positive, top_negative)
    
scatter1 <- ggplot(data = vol_val, aes(x = scaled, y = -log10(FDR))) +
  geom_point(aes(color = color), size = 3.5) +
  geom_text_repel(aes(label = motif), size = 4, max.overlaps = 20, 
                  data = top_20) +  # Label only the selected 20 points
  labs(title = paste0(celltype),
       x = "Age correlated activity (scaled pcc)",
       y = "-Log10 FDR") +
  theme_classic() +
  theme(text = element_text(size = 20)) +
#  geom_hline(yintercept = -log10(0.05), linetype = "longdash", color = "black") +
  geom_vline(xintercept = 0.0, linetype = "longdash", color = "black") +
#  geom_vline(xintercept = -0.05, linetype = "longdash", color = "black") +
  coord_cartesian(expand = TRUE) +
#  scale_x_continuous(limits = c(-3, 3)) +
#  scale_y_continuous(limits = c(0, )) +
  scale_color_identity()


assign(paste0('plot_', celltype), scatter1)
assign(paste0('vol_val_', celltype), vol_val)

write.table(vol_val, file = paste0('motif_atac_pcc/', celltype, '_all_motifs_fdr_scaled_pcc_prox_dist_ratio.tsv'), sep = '\t', quote=F)

    }

plot_Oligo 
plot_OPC 
plot_SUB 
plot_CA1 
plot_Macro 
plot_Endo 
plot_VLMC 
`plot_CA2-CA3`
plot_SST 
plot_VIP 
plot_LAMP5 
plot_Chandelier 
plot_NR2F2 
plot_PVALB 
`plot_T-Cell` 
plot_Microglia 
plot_Astro 
plot_DG
