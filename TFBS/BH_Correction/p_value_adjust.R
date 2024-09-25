# Load necessary libraries
library(dplyr)

# Ensure data is loaded (replace 'your_file.tsv' with your actual file path)
data_filtered <- read.delim('contingency_fimo_cub.tsv', header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Ensure 'Gene' column is ignored and filter out any unwanted columns
data_filtered <- data_filtered %>% select(-Gene)

# Convert TFBS columns, CAI, and Nc to numeric (if not already)
tfbs_columns <- colnames(data_filtered)[!colnames(data_filtered) %in% c("CAI", "Nc")]
data_filtered[tfbs_columns] <- lapply(data_filtered[tfbs_columns], function(x) as.numeric(as.character(x)))
data_filtered$CAI <- as.numeric(as.character(data_filtered$CAI))
data_filtered$Nc <- as.numeric(as.character(data_filtered$Nc))

# Remove columns with zero variance
zero_variance_tfbs <- sapply(data_filtered[tfbs_columns], function(x) var(x, na.rm = TRUE) == 0)
if (any(zero_variance_tfbs)) {
  data_filtered <- data_filtered[, !zero_variance_tfbs]
  tfbs_columns <- colnames(data_filtered)[!colnames(data_filtered) %in% c("CAI", "Nc")]
}

# Initialize dataframes to store correlation results for CAI and Nc separately
results_cai <- data.frame(TFBS = character(), Correlation = numeric(), P_value = numeric(), stringsAsFactors = FALSE)
results_nc <- data.frame(TFBS = character(), Correlation = numeric(), P_value = numeric(), stringsAsFactors = FALSE)

# Perform correlation for each TFBS with CAI and Nc using Spearman correlation
for (tfbs in tfbs_columns) {
  # Correlation with CAI
  corr_cai <- tryCatch({
    cor.test(data_filtered[[tfbs]], data_filtered$CAI, method = "spearman")
  }, error = function(e) {
    message(paste("Error in correlation between", tfbs, "and CAI:", e$message))
    return(NULL)
  })
  
  if (!is.null(corr_cai)) {
    results_cai <- rbind(results_cai, data.frame(TFBS = tfbs, Correlation = corr_cai$estimate, P_value = corr_cai$p.value))
  }
  
  # Correlation with Nc
  corr_nc <- tryCatch({
    cor.test(data_filtered[[tfbs]], data_filtered$Nc, method = "spearman")
  }, error = function(e) {
    message(paste("Error in correlation between", tfbs, "and Nc:", e$message))
    return(NULL)
  })
  
  if (!is.null(corr_nc)) {
    results_nc <- rbind(results_nc, data.frame(TFBS = tfbs, Correlation = corr_nc$estimate, P_value = corr_nc$p.value))
  }
}

# Apply Benjamini-Hochberg correction to CAI p-values
results_cai$Adjusted_P_value <- p.adjust(results_cai$P_value, method = "holm")

# Add a new column to reject or accept the null hypothesis based on the adjusted p-value
results_cai$Null_Hypothesis <- ifelse(results_cai$Adjusted_P_value < 0.05, "Reject", "Accept")

# Apply Benjamini-Hochberg correction to Nc p-values
results_nc$Adjusted_P_value <- p.adjust(results_nc$P_value, method = "holm")
results_nc$Null_Hypothesis <- ifelse(results_nc$Adjusted_P_value < 0.05, "Reject", "Accept")

# Save the full CAI results to a CSV file
write.csv(results_cai, "correlation_results_cai_all.csv", row.names = FALSE)

# Save the full Nc results to a CSV file
write.csv(results_nc, "correlation_results_nc_all.csv", row.names = FALSE)

# Filter for statistically significant CAI p-values (< 0.05)
results_cai_significant <- results_cai %>% filter(P_value < 0.05)

# Apply holm correction again only on the significant CAI p-values
results_cai_significant$Adjusted_P_value <- p.adjust(results_cai_significant$P_value, method = "holm")
results_cai_significant$Null_Hypothesis <- ifelse(results_cai_significant$Adjusted_P_value < 0.05, "Reject", "Accept")

# Save the significant CAI results to a CSV file
write.csv(results_cai_significant, "correlation_results_cai_significant.csv", row.names = FALSE)

# Filter for statistically significant Nc p-values (< 0.05)
results_nc_significant <- results_nc %>% filter(P_value < 0.05)

# Apply holm correction again only on the significant Nc p-values
results_nc_significant$Adjusted_P_value <- p.adjust(results_nc_significant$P_value, method = "holm")
results_nc_significant$Null_Hypothesis <- ifelse(results_nc_significant$Adjusted_P_value < 0.05, "Reject", "Accept")

# Save the significant Nc results to a CSV file
write.csv(results_nc_significant, "correlation_results_nc_significant.csv", row.names = FALSE)

# Read the CSV files into data frames
df1 <- read.csv("correlation_results_cai_significant.csv")
df2 <- read.csv("correlation_results_nc_significant.csv")

# Merge the data frames based on the 'TFBS' column
merged_df <- merge(df1, df2, by = "TFBS")

# Rename the columns to append '_CAI' and '_Nc' instead of '.x' and '.y'
colnames(merged_df) <- gsub("\\.x$", "_CAI", colnames(merged_df))
colnames(merged_df) <- gsub("\\.y$", "_Nc", colnames(merged_df))

# Write the merged data frame to a new CSV file
write.csv(merged_df, "common_significant_TFBS.csv", row.names = FALSE)

