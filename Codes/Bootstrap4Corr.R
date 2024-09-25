#WRKY CUI analysis

#importing necessary libraries
library(tidyverse)
library(corrplot)
library(readxl)
library(writexl)
library(ggpmisc)
library(boot)

#loading the data
df_plot <- read_excel("cub_omega.xlsx")

#getting subset of the data
omega_less <- subset(df_plot, w_ratio >1 & w_ratio < 400)

# Compute Spearman correlation
result <- cor.test(omega_less$w_ratio, omega_less$CAI, method = "spearman")

# Print results
print(result$estimate)  # Spearman correlation coefficient
print(result$p.value)   # p-value


#correlation matrix for purifying selection 
corrmtx_1 <- cor(omega_less[-c(1,4,5)],method=c('spearman'))
write.table(corrmtx_1,
            file = 'corr_coef_1.txt',
            sep = '\t')
p_val <- cor.mtest(omega_less[-c(1,4,5)],conf.level=0.95)
write.table(p_val,
            file = 'corr_coef_1.txt',
            sep = '\t',
            append = TRUE) 
#
#Correlation matrix plot
corrplot(corrmtx_1, 
         p.mat = p_val$p, 
         method = 'color', 
         sig.level = c(0.001, 0.01, 0.05),
         diag = FALSE,
         insig = 'label_sig', 
         order = 'hclust',
         addrect = 2,
         pch.cex = 1,
         tl.col = 'black')

#Neutrality plot
ggplot(omega_less, aes(CAI, w_ratio)) +
  geom_point() +
  geom_smooth(method = "lm")+
  stat_poly_eq(use_label(c('eq','R2','p.value')))+
  theme_classic()

#install.packages("boot")

new_dataframe <- data.frame(w_ratio = omega_less$w_ratio, CAI = omega_less$CAI)

# Assuming you have created 'new_dataframe' with two columns: 'w_ratio' and 'CAI'

# Define the correlation function for bootstrapping
correlation_fn <- function(data, indices) {
  # Sample the data based on the bootstrap indices
  d <- data[indices, ]
  # Return the Pearson correlation between 'w_ratio' and 'CAI'
  return(cor(d$w_ratio, d$CAI, method = "spearman"))  # You can also use method = "spearman" if needed
}

# Perform the bootstrapping to estimate the correlation coefficient
bootstrap_results <- boot(data = new_dataframe, statistic = correlation_fn, R = 1000)

# View the bootstrapping results
print(bootstrap_results)

# Compute the confidence intervals
ci <- boot.ci(bootstrap_results, type = "perc")  # Use "bca" or "norm" as alternative methods
print(ci)


