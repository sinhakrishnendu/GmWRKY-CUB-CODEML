library(readxl)
library(Hmisc)
library(ggplot2)
library(reshape2)
library(gridExtra)

# Load the Excel file and specify the sheet name
file_path <- "cub_omega.xlsx"  # replace with your file path
sheet_name <- "Sheet9"  # replace with your sheet name
data <- read_excel(file_path, sheet = sheet_name)

# Keep only numeric columns
numeric_data <- data[sapply(data, is.numeric)]

# Calculate the correlation matrix and p-values
correlation_results <- rcorr(as.matrix(numeric_data), type = "spearman")

# Extract correlation matrix and p-values
correlation_matrix <- correlation_results$r
p_value_matrix <- correlation_results$P

# Combine correlation coefficients and p-values into a single matrix
combined_matrix <- correlation_matrix
combined_matrix[] <- sprintf("%.2f\n(%.2f)", correlation_matrix, p_value_matrix)

# Convert the combined matrix to a data frame for ggplot
combined_df <- melt(as.matrix(combined_matrix), varnames = c("X1", "X2"), value.name = "Value")
combined_df$X1 <- factor(combined_df$X1, levels = rev(colnames(combined_matrix)))
combined_df$X2 <- factor(combined_df$X2, levels = rev(colnames(combined_matrix)))

# Plot the combined matrix using ggplot2
p <- ggplot(combined_df, aes(X2, X1, fill = as.numeric(gsub(".*\\(", "", Value)))) +
  geom_tile() +
  geom_text(aes(label = Value), size = 4) +
  scale_fill_gradient2(low = "red", high = "blue", mid = "white", midpoint = 0) +
  labs(x = "", y = "", fill = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display the plot
print(p)
