# Load necessary library
install.packages('openxlsx')
library(openxlsx)

# Read the TSV file
data <- read.table("fimo.tsv", header = TRUE, sep = "\t")

# Create a contingency table
contingency_table <- table(data$sequence_name, data$motif_alt_id)

# Convert the contingency table to a matrix
pivot_matrix <- as.matrix(contingency_table)

colnames(pivot_matrix)[0] <- "Name"

# Write the matrix to an Excel file
write.xlsx(pivot_matrix, file = "TF_GENE_matrix.xlsx")

# Load necessary libraries
library(readxl)   # For reading Excel files
library(pheatmap) # For creating heatmaps

# Read the Excel file
# Replace 'your_file.xlsx' with the path to your Excel file and 'Sheet1' with the sheet name if applicable
data <- read_excel('TF_GENE_matrix.xlsx')

# Extract the column names for x-axis
x_axis_labels <- data[[1]]

# Remove the first column (Names) and convert the remaining data to matrix
data_matrix <- as.matrix(data[,-1])

# Set row names to the first column values
rownames(data_matrix) <- x_axis_labels

# Create a heatmap
pheatmap(data_matrix, 
         cluster_rows = TRUE,   # Disable clustering rows
         cluster_cols = TRUE,   # Disable clustering columns
         main = 'Heatmap of Contingency Table',
         scale = 'none',
         fontsize = 4,    # Adjust font size for better readability
         fontsize_row = 4, 
         fontsize_col = 4)         # No scaling
