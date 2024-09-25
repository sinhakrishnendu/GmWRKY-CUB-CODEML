#WRKY CUI analysis

#importing necessary libraries
library(tidyverse)
library(corrplot)
library(readxl)
library(writexl)
library(ggpmisc)

#reading codon usage dataset from codonw
df_plot <- read_excel("cub_omega1.xlsx")
#View(cui_dataset)

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
#
#Neutrality plot
ggplot(df_plot, aes(CAI, w_ratio)) +
  geom_point() +
  geom_smooth(method = "lm")+
  stat_poly_eq(use_label(c('eq','R2','p.value')))+
  theme_classic()
#
#regression data of neutrality plot
model <- lm(CAI~w_ratio,data = df_plot)
summary(model)
#
GC3s <- seq(0,1,length.out = 10000)
Nc_e <- 2+GC3s+(29/(GC3s^2+(1-GC3s)^2))
expEnc <- data.frame(GC3s,Nc_e)
#
#Nc plot
ggplot()+
  geom_point(data = df_plot, aes(x = GC3s, y = Nc)) +
  geom_line(data = expEnc,aes(x = GC3s,y = Nc_e))+
  theme_classic()
#
library(ggrepel)
#PR2 plot
ggplot()+
  geom_point(data = df_plot,aes(x=PR2_x,y=PR2_y))+
  geom_vhlines(xintercept = 0.5, yintercept = 0.5)+
  xlim(0,1) + ylim(0,1)+
  xlab('G3/(G3+C3)') + ylab('A3/(A3+T3)')+
  theme_classic()
#
#retrieving mean value of dependent and independent variables
summary(df_plot$PR2_x)
summary(df_plot$PR2_y)
#
#hist
ggplot()+
  geom_histogram(data = df_plot,aes(x = eNc_Nc_hist),
                 binwidth = 0.02,color = 'white')+
  xlab('(Nc_exp-Nc_obs)/Nc_exp')+
  ylab('Frequency')+
  theme_classic()
#
hist_data <- hist(df_plot$eNc_Nc_hist, 
                  xlab = '(eNc-Nc)/eNc',
                  main = '', breaks = 5,
                  ylim = c(0,70),
                  border = 'white',
                  col = 'black')
text(hist_data$mids, hist_data$counts, 
     labels = round(hist_data$counts / sum(hist_data$counts) * 100,1), 
     adj = c(0.5, -0.5))
#
#Kolmogorov-Smirnov test
ks.test(df_plot$Nc,
        expEnc$Nc_e,
        alternative = 'less')


df_plot_1 <- read_excel("cub_omega.xlsx")

omega_less_1 <- subset(df_plot_1, w_ratio >1 & w_ratio < 4)
#correlation matrix for purifying selection 
corrmtx_1 <- cor(omega_less_1[-c(1,4,5)],method=c('spearman'))
write.table(corrmtx_1,
            file = 'corr_coef_1.txt',
            sep = '\t')
p_val <- cor.mtest(omega_less_1[-c(1,4,5)],conf.level=0.95)
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
ggplot(omega_less_1, aes(CAI, w_ratio)) +
  geom_point() +
  geom_smooth(method = "lm")+
  stat_poly_eq(use_label(c('eq','R2','p.value')))+
  theme_classic()

#outlier detection
#install.packages("outliers")
#library(outliers)

# Detect outliers using Grubbs' test
#grubbs.test(omega_less_1$w_ratio)

#boxplot(omega_less_1$w_ratio, main = "Boxplot for Outlier Detection")

# Identify outliers using IQR method
#Q1 <- quantile(omega_less_1$w_ratio, 0.25)
#Q3 <- quantile(omega_less_1$w_ratio, 0.75)
#IQR <- Q3 - Q1

# Thresholds for detecting outliers
#lower_bound <- Q1 - 1.5 * IQR
#upper_bound <- Q3 + 1.5 * IQR

# Outliers
#outliers <- omega_less_1$w_ratio[omega_less_1$w_ratio < lower_bound | omega_less_1$w_ratio > upper_bound]
#outliers
# Install and load the boot package
#install.packages("boot")
library(boot)
new_dataframe <- data.frame(w_ratio = omega_less_1$w_ratio, CAI = omega_less_1$CAI)

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


