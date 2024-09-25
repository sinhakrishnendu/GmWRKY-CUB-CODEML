# Load required libraries
library(ggplot2)
library(reshape2)
library(ggcorrplot)

# Load the data
data <- read.csv("tfbs_cub_all_significant_output.tsv", sep = "\t", header = TRUE)

# Drop the 'Gene' column (since it's not numeric and can't be correlated)
data_numeric <- data[, -which(names(data) == "Gene")]

# Select only 'CAI', 'Nc', and TF columns
cai_nc_columns <- data_numeric[, c('CAI', 'Nc')]
tf_columns <- data_numeric[, !(names(data_numeric) %in% c('CAI', 'Nc'))]

# Compute the correlation matrix for TFs against CAI and Nc only
corr_matrix <- cor(cbind(cai_nc_columns, tf_columns), method = "spearman", use = "complete.obs")

# Subset the correlation matrix to only include correlations between CAI/Nc and TFs
corr_matrix_subset <- corr_matrix[c('CAI', 'Nc'), !(colnames(corr_matrix) %in% c('CAI', 'Nc'))]

# Convert the subsetted correlation matrix to long format for plotting
corr_data_long <- melt(corr_matrix_subset)

# Create a correlation heatmap plot without R values and setting the color range
ggplot(corr_data_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", high = "blue", mid = "white", 
                       midpoint = 0, limit = c(-0.4, 0.4), 
                       name = "Spearman\nCorrelation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),  # Align x-axis labels properly
        axis.text.y = element_text(size = 8)) +              # Adjust y-axis label size
  labs(title = "Correlation Between Transcription Factors and CAI/Nc",
       x = "CAI/Nc", y = "Transcription Factors")

# Save the plot to a file
ggsave("tf_cai_nc_correlation_plot.png", width = 10, height = 8)
#
#
#heatmap
# Load required libraries
# Load required libraries
# Load required libraries
# Load required libraries
library(tidyverse) # For data manipulation
library(pheatmap)  # For heatmap plotting

# Define the significance threshold
threshold <- 0.05  # Adjust this value as needed for your analysis

# Step 1: Read the data from the TSV file
data <- read_tsv("tfbs_cub_all_significant_output.tsv")

# Step 2: Drop the 'Gene' column
data_clean <- data %>%
  select(-Gene)

# Separate TFBS frequencies, CAI, and Nc
tfbs_data <- data_clean %>%
  select(-CAI, -Nc)

ca_nc_data <- data_clean %>%
  select(CAI, Nc)

# Step 3: Compute the Spearman correlation between TFBS frequencies and CAI/Nc
# Calculate Spearman correlations
correlation_matrix <- cor(tfbs_data, ca_nc_data, method = "spearman")

# Step 4: Replace insignificant values with NA
correlation_matrix[abs(correlation_matrix) < threshold] <- NA

# Step 5: Remove rows and columns with too many NA values
na_threshold <- 0.5 # Define the maximum proportion of NA values allowed
na_counts <- apply(is.na(correlation_matrix), 1, mean)
valid_rows <- na_counts < na_threshold

na_counts <- apply(is.na(correlation_matrix), 2, mean)
valid_cols <- na_counts < na_threshold

cleaned_matrix <- correlation_matrix[valid_rows, valid_cols]

# Step 6: Create a high-resolution heatmap with hierarchical clustering
output_file <- "heatmap_correlation_high_res.png"

png(output_file, width = 1200, height = 1600, res = 300) # High resolution
pheatmap(cleaned_matrix,
         main = "Spearman Correlation Heatmap",
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = FALSE, # Omit R values
         fontsize_row = 6, # Small font size for TFBS names
         fontsize_col = 10, # Font size for CAI and Nc labels
         angle_row = 45,   # Tilt only TFBS names
         na_col = "orange")  # Color for NA (blank) cells
dev.off() # Close the PNG device



