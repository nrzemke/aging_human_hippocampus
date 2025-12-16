# Extract expression values and metadata
cdkn2a_expr <- GetAssayData(Human_Combined_Log_Sctransformed_Subset, assay = "SCT", slot = "data")["CDKN2A", ]
meta <- Human_Combined_Log_Sctransformed_Subset@meta.data

# Ensure correct order
stopifnot(all(names(cdkn2a_expr) == rownames(meta)))

# Add expression to metadata
meta$CDKN2A_expr <- cdkn2a_expr

# Split by Final1 group
meta_split <- split(meta, meta$Final1)

# Load age info
age <- read.table('Age_donor_ids.tsv', header = TRUE)

# Load library for join
library(plyr)  # or dplyr if you prefer

# Function to compute percentage
calc_percent <- function(df, celltype) {
  out <- aggregate(CDKN2A_expr >= 1 ~ orig.ident, data = df, FUN = function(x) mean(x) * 100)
  colnames(out)[2] <- "Percent_CDKN2A_Expr_1plus"
  out$celltype <- celltype
  return(out)
}

# Loop to calculate, join age, and write file
for (final1_group in names(meta_split)) {
  df_group <- meta_split[[final1_group]]
  
  # Calculate CDKN2A expression %
  result <- calc_percent(df_group, final1_group)
  
  # Join age info
  result <- join(result, age, by = 'orig.ident')  # If using dplyr: left_join(result, age, by = 'orig.ident')
  
  # Create safe filename
  file_name <- paste0("CDKN2A_Percent_", gsub("[^A-Za-z0-9_]", "_", final1_group), ".tsv")
  
  # Write to file
  write.table(result, file = file_name, sep = "\t", row.names = FALSE, quote = FALSE)
}
