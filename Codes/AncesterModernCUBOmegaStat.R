#WRKY CUI analysis

#importing necessary libraries
library(tidyverse)
library(corrplot)
library(readxl)
library(writexl)
library(ggpmisc)
library(ggpubr)

#reading codon usage dataset from codonw
df_ancestral <- read_excel("AncestralCUBOmega.xlsx")
df_extant <- read_excel('ExtantCUBOmega.xlsx')

###
#subsetting dataframe
df_asr <- subset(df_ancestral, w_ratio < 10)
df_wrky <- subset(df_extant, w_ratio < 10)

#checking normality assumption on ASR CUB
shapiro.test(df_ancestral$CAI)
shapiro.test(df_asr$CAI)
shapiro.test(df_ancestral$Nc)
shapiro.test(df_asr$Nc)

# Calculate Spearman correlation and p-value
cor_test <- cor.test(df_asr$CAI, df_asr$Nc, method = "spearman")
spearman_corr <- cor_test$estimate
p_value <- cor_test$p.value

# Print Spearman correlation and p-value
print(paste("Spearman correlation:", round(spearman_corr, 3)))
print(paste("p-value:", format(p_value, scientific = TRUE)))

#regression plot
ggplot(df_asr, aes(CAI, Nc)) +
  geom_point() +
  geom_smooth(method = "lm")+
  stat_poly_eq(use_label(c('eq','R2','p.value')))+
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 16),  # axis lable size
    axis.title.y = element_text(size = 16)   # do
  )
ggsave('RegPlot.tiff',dpi = 600)


###

###
#Mann-Whitney U
#for CAI
mwu_cai <- wilcox.test(df_asr$CAI, df_wrky$CAI) 
print(mwu_cai)
#for Nc
mwu_nc <- wilcox.test(df_asr$Nc, df_wrky$Nc)
print(mwu_nc)
#for omega ratio
mwu_omega <- wilcox.test(df_asr$w_ratio,df_wrky$w_ratio)
print(mwu_omega)
###

###
#Graphical representation
# CAI
df_combined <- data.frame(
  CAI = c(df_asr$CAI, df_wrky$CAI),
  Group = factor(c(rep("Ancestral", length(df_asr$CAI)), rep("Modern", length(df_wrky$CAI))))
)
cai_v <- ggviolin(df_combined, x = 'Group', y = 'CAI', fill = 'Group',
         palette = c("#00AFBB", "#E7B800"),
         add = "boxplot", add.params = list(fill = "#f5f4f0"))+
         labs(x = "Gene Type")+
         stat_compare_means(comparisons = list(c("Ancestral", "Modern")),
                            method = "wilcox.test", 
                            label = "p.signif")
cai_v

ggsave("violin_plot_mann_whitney_CAI.tiff", plot = cai_v, width = 4, height = 4, dpi = 900)# Save the plot as a high-resolution image

#Nc
df_combined_nc <- data.frame(
  Nc = c(df_asr$Nc, df_wrky$Nc),
  Group = factor(c(rep("Ancestral", length(df_asr$Nc)), rep("Modern", length(df_wrky$Nc))))
)

nc_v <- ggviolin(df_combined_nc, x = 'Group', y = 'Nc', fill = 'Group',
                  palette = c("#00AFBB", "#E7B800"),
                  add = "boxplot", add.params = list(fill = "#f5f4f0"))+
  labs(x = "Gene Type")+
  stat_compare_means(comparisons = list(c("Ancestral", "Modern")),
                     method = "wilcox.test", 
                     label = "p.signif")
nc_v

ggsave("violin_plot_mann_whitney_Nc.tiff", plot = nc_v, width = 4, height = 4, dpi = 900)# Save the plot as a high-resolution image

#omega
df_combined_w_ratio <- data.frame(
  w_ratio = c(df_asr$w_ratio, df_wrky$w_ratio),
  Group = factor(c(rep("Ancestral", length(df_asr$w_ratio)), rep("Modern", length(df_wrky$w_ratio))))
)

w_v <- ggviolin(df_combined_w_ratio, x = 'Group', y = 'w_ratio', fill = 'Group',
                   palette = c("#00AFBB", "#E7B800"),
                   add = "boxplot", add.params = list(fill = "#f5f4f0"))+
  labs(x = "Gene Type", y = expression(omega)) + 
  stat_compare_means(comparisons = list(c("Ancestral", "Modern")),
                     method = "wilcox.test", 
                     label = "p.signif")
w_v

ggsave("violin_plot_mann_whitney_w.tiff", plot = w_v, width = 4, height = 4, dpi = 900)# Save the plot as a high-resolution image
###
















df_combined <- data.frame(
  CAI = c(df_asr$CAI, df_wrky$CAI),
  Group = factor(c(rep("Ancestral", length(df_asr$CAI)), rep("Modern", length(df_wrky$CAI))))
)
cai_v <- ggviolin(df_combined, x = 'Group', y = 'CAI', fill = 'Group',
                  palette = c("#aba89f", "#a12d2d"),
                  add = "boxplot", add.params = list(fill = "#f5f4f0"))+
  labs(x = "Gene Type")+
  stat_compare_means(comparisons = list(c("Ancestral", "Modern")),
                     method = "wilcox.test", 
                     label = "p.signif")
cai_v
ggsave('x.tiff',plot = cai_v, width = 4, height = 4, dpi = 900)
