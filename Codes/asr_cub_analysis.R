#WRKY CUI analysis

#importing necessary libraries
library(tidyverse)
library(corrplot)
library(readxl)
library(writexl)
library(ggpmisc)

#reading codon usage dataset from codonw
df_asr <- read_excel("ASR_CUB.xlsx")
df_gmwrky <- read_excel('CUBindices.xlsx')
df_omega <- read_excel('cub_omega.xlsx')

###
# Calculate Spearman correlation and p-value
cor_test <- cor.test(df_plot$CAI, df_plot$Nc, method = "spearman")
spearman_corr <- cor_test$estimate
p_value <- cor_test$p.value

# Print Spearman correlation and p-value
print(paste("Spearman correlation:", round(spearman_corr, 3)))
print(paste("p-value:", format(p_value, scientific = TRUE)))

#linear model
ggplot(df_plot, aes(CAI, Nc)) +
  geom_point() +
  geom_smooth(method = "lm")+
  stat_poly_eq(use_label(c('eq','R2','p.value')))+
  theme_classic()


#correlation matrix
corrmtx <- cor(df_plot[-c(1)],method=c('spearman'))
write.table(corrmtx,
            file = 'corr_coef.txt',
            sep = '\t')
p_val <- cor.mtest(df_plot[-c(1)],conf.level=0.95)
write.table(p_val,
            file = 'corr_coef.txt',
            sep = '\t',
            append = TRUE) 
#
#Correlation matrix plot
corrplot(corrmtx, 
         p.mat = p_val$p, 
         method = 'color', 
         sig.level = c(0.001, 0.01, 0.05),
         diag = FALSE,
         insig = 'label_sig', 
         order = 'hclust',
         addrect = 2,
         pch.cex = 1,
         tl.col = 'black')

#Mann-Whitney U
mann_whitney <- wilcox.test(df_asr$CAI,df_gmwrky$CAI)
print(mann_whitney)

mann_whitney <- wilcox.test(df_asr$Nc,df_gmwrky$Nc)
print(mann_whitney)

median(df_asr$Nc)
median(df_gmwrky$Nc)

#
median(df_asr$CAI)
median(df_gmwrky$CAI)

# Install and load the necessary package
install.packages("dunn.test")
library(dunn.test)

# Combine the data into a single data frame
data <- data.frame(
  value = c(df_asr$CAI, df_asr$Nc, df_gmwrky$CAI, df_gmwrky$Nc),
  group = factor(c(rep("df_asr_CA", length(df_asr$CAI)),
                   rep("df_asr_Nc", length(df_asr$Nc)),
                   rep("df_gmwrky_CA", length(df_gmwrky$CAI)),
                   rep("df_gmwrky_Nc", length(df_gmwrky$Nc))))
)

# Perform Kruskal-Wallis test
kruskal_result <- kruskal.test(value ~ group, data = data)
print(kruskal_result)

# Perform Dunn's test for post-hoc analysis
dunn_result <- dunn.test(data$value, data$group, method = "bonferroni")
print(dunn_result)

#t-test
t <- t.test(df_asr$CAI,df_gmwrky$CAI)
print(t)

###
###

#Mann-Whitney U
mann_whitney <- wilcox.test(subset(df_asr, w_ratio <10)$w_ratio,subset(df_omega,w_ratio < 10)$w_ratio)
print(mann_whitney)

mann_whitney <- wilcox.test(df_asr$Nc,df_gmwrky$Nc)
print(mann_whitney)

median(subset(df_asr, w_ratio <10)$w_ratio)
median(subset(df_omega,w_ratio < 10)$w_ratio)
mean(subset(df_asr, w_ratio <10)$w_ratio)
mean(subset(df_omega,w_ratio < 10)$w_ratio)
#
median(df_asr$CAI)
median(df_gmwrky$CAI)

# Install and load the necessary package
install.packages("dunn.test")
library(dunn.test)

# Combine the data into a single data frame
data <- data.frame(
  value = c(df_asr$CAI, df_asr$Nc, df_gmwrky$CAI, df_gmwrky$Nc),
  group = factor(c(rep("df_asr_CA", length(df_asr$CAI)),
                   rep("df_asr_Nc", length(df_asr$Nc)),
                   rep("df_gmwrky_CA", length(df_gmwrky$CAI)),
                   rep("df_gmwrky_Nc", length(df_gmwrky$Nc))))
)

# Perform Kruskal-Wallis test
kruskal_result <- kruskal.test(value ~ group, data = data)
print(kruskal_result)

# Perform Dunn's test for post-hoc analysis
dunn_result <- dunn.test(data$value, data$group, method = "bonferroni")
print(dunn_result)

#t-test
t <- t.test(subset(df_asr, w_ratio <10)$w_ratio,subset(df_omega,w_ratio < 10)$w_ratio)
print(t)

#x<-subset(df_asr, w_ratio <100 & w_ratio >1)

shapiro.test((subset(df_asr, w_ratio <10))$w_ratio)
shapiro.test((subset(df_omega, w_ratio <10))$w_ratio)

ks.test((subset(df_asr, w_ratio <10))$w_ratio, "pnorm", mean((subset(df_asr, w_ratio <10))$w_ratio), sd((subset(df_asr, w_ratio <10))$w_ratio))
ks.test(subset(df_asr, w_ratio <10)$w_ratio,subset(df_omega,w_ratio < 10)$w_ratio,alternative = 'greater')

#plotting cdfs
# Subset the data for w_ratio < 10
df_asr_subset <- subset(df_asr, w_ratio < 10)
df_omega_subset <- subset(df_omega, w_ratio < 10)

# Combine the data into a single data frame
df_combined <- rbind(
  data.frame(w_ratio = df_asr_subset$w_ratio, group = "df_asr"),
  data.frame(w_ratio = df_omega_subset$w_ratio, group = "df_omega")
)

# Plot the ECDFs using ggplot2
ggplot(df_combined, aes(x = w_ratio, color = group)) +
  stat_ecdf(geom = "step", size = 1) +
  labs(
    title = "Empirical Cumulative Distribution Functions (ECDF)",
    x = "w_ratio",
    y = "ECDF"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank()
  )

