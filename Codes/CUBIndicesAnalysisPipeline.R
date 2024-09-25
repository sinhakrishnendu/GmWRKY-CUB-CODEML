#WRKY CUI analysis

#importing necessary libraries
library(tidyverse)
library(corrplot)
library(readxl)
library(writexl)
library(ggpmisc)
library(ggrepel)

#reading codon usage dataset from codonw
cui_dataset <- read_excel("CUI.xlsx")
#View(cui_dataset)

#reading codon usage dataset from CAIcal
caical <- read_excel("CAICal_GC123.xlsx")
#View(caical)

#merging and cleaning two datasets to get final dataset form further analysis based on TF IDs

joined_df <- merge(cui_dataset,caical[c(1,25,39,53,49:52)],
                   by.x = 'title',by.y = 'sequences',
                   all.x = TRUE,all.y = TRUE)
view(joined_df)

#creating new columns with unique values 
joined_df$GC1 <- joined_df$`%G1+C1`/100
joined_df$GC2 <- joined_df$`%G2+C2`/100
joined_df$GC3 <- joined_df$`%G3+C3`/100
joined_df$T3 <- joined_df$`%T3`/100
joined_df$C3 <- joined_df$`%C3`/100
joined_df$A3 <- joined_df$`%A3`/100
joined_df$G3 <- joined_df$`%G3`/100
view(joined_df)

#generating final data frame for analysis
cui <- joined_df[-c(7:8,12:22)]
View(cui)
write.csv(cui,file = 'Supplimentary_datafile.csv',row.names = FALSE)

#checking normality of the variables by Shapiro-Wilk normality test

shapiro_wilk <- function(x){
  shapiro.test(x)$p.value} #creating function to perform Shapiro-Wilk test and return the p-values only
#
shapiro_result <- apply(cui[-c(1)],2,shapiro_wilk) #applying the function to every numerical column(2)
#
write.csv(shapiro_result,file='shapiro_test_result.csv',
          row.names = colnames(cui[-c(1)]))#exporting the output
#
print(shapiro_result)
#
#dataframe for plots
df_plot <- data.frame(title=cui$title)
df_plot$GC12 <- (cui$GC1+cui$GC2)/2
df_plot$GC3 <- cui$GC3
df_plot$Nc <- cui$Nc
df_plot$Nc_exp <- 2+cui$GC3s+(29/(cui$GC3s^2+(1-cui$GC3s)^2))
df_plot$GC3s <- cui$GC3s
df_plot$PR2_y <- cui$A3/(cui$A3+cui$T3)
df_plot$PR2_x <- cui$G3/(cui$G3+cui$C3)
df_plot$eNc_Nc_hist <- (df_plot$Nc_exp-df_plot$Nc)/df_plot$Nc_exp
#
View(df_plot)
write.csv(df_plot,file = 'df_for_plots.csv')
#
#correlation matrix
corrmtx <- cor(cui[-c(1)],method=c('spearman'))
write.table(corrmtx,
            file = 'corr_coef.txt',
            sep = '\t')
p_val <- cor.mtest(cui[-c(1)],conf.level=0.95)
write.table(p_val,
            file = 'corr_coef.txt',
            sep = '\t',
            append = TRUE) 
###
tiff("correlation_plot.tiff", width = 1920, height = 1080, res = 300)
#Correlation matrix plot
corrplot(corrmtx, 
         p.mat = p_val$p, 
         method = 'color', 
         sig.level = c(0.001, 0.01, 0.05),
         diag = FALSE,
         insig = 'label_sig', 
         order = 'hclust',
         addrect = 2,
         pch.cex = 0.7,
         tl.col = 'black',
         col = COL2('RdBu', 20))
dev.off()
###

###
#Neutrality plot
ggplot(df_plot, aes(GC3, GC12)) +
  geom_point() +
  geom_smooth(method = "lm")+
  stat_poly_eq(use_label(c('eq','R2','p.value')))+
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 16),  # axis lable size
    axis.title.y = element_text(size = 16)   # do
  )
ggsave('NeutralityPlot.tiff',dpi = 600)

#regression data of neutrality plot
model <- lm(GC12~GC3,data = df_plot)
summary(model)
###

###
#expected ENC
GC3s <- seq(0,1,length.out = 179)
Nc_e <- 2+GC3s+(29/(GC3s^2+(1-GC3s)^2))
expEnc <- data.frame(GC3s,Nc_e)

#Nc plot
ggplot()+
  geom_point(data = df_plot, aes(x = GC3s, y = Nc)) +
  geom_line(data = expEnc,aes(x = GC3s,y = Nc_e))+
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 16),  # axis lable size
    axis.title.y = element_text(size = 16)   # do
  )

ggsave('NcPlot.tiff',dpi = 300)
###

###
#hist
ggplot()+
  geom_histogram(data = df_plot,aes(x = eNc_Nc_hist),
                 binwidth = 0.02,color = 'white')+
  xlab('(Nc_exp-Nc_obs)/Nc_exp')+
  ylab('Frequency')+
  theme_classic()

tiff("eNcHist.tiff", width = 1080, height = 1920, res = 300)

hist_data <- hist(df_plot$eNc_Nc_hist, 
                  xlab = '(eNc-Nc)/eNc',
                  main = '', breaks = 5,
                  ylim = c(0,70),
                  border = 'white',
                  col = 'black')
text(hist_data$mids, hist_data$counts, 
     labels = round(hist_data$counts / sum(hist_data$counts) * 100,1), 
     adj = c(0.5, -0.5))
dev.off()
###

###
# PR2 plot
ggplot() +
  geom_point(data = df_plot, aes(x = PR2_x, y = PR2_y)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "darkred") +  # Corrected geom_hline
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "darkred") +  # Corrected geom_vline
  xlim(0, 1) + ylim(0, 1) +
  xlab('G3/(G3+C3)') + ylab('A3/(A3+T3)') +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16),  # axis lable size
    axis.title.y = element_text(size = 16)   # do
  )

# Save the plot
ggsave('PR2Plot.tiff', dpi = 300)

#retrieving mean value of dependent and independent variables
summary(df_plot$PR2_x)
summary(df_plot$PR2_y)
###

###
#Kolmogorov-Smirnov test
ks.test(df_plot$Nc,
        expEnc$Nc_e,
        alternative = 'less')
###

###
#KS CDF visualization
data_plot <- data.frame(
  value = c(df_plot$Nc, expEnc$Nc_e),
  group = rep(c("Nc", "eNc"), times = c(length(df_plot$Nc), length(expEnc$Nc_e)))
)

# Create the CDF plot
cdf_ks <- ggplot(data_plot, aes(x = value, color = group)) +
  stat_ecdf(geom = "step") +
  labs(title = "Cumulative Distribution Functions",
       x = "Value",
       y = "CDF") +
  scale_color_manual(values = c("Nc" = "blue", "eNc" = "red")) +
  theme(
    panel.background = element_rect(fill = "white"), # White background
    panel.grid = element_blank(),                    # No gridlines
    axis.line = element_line(color = "black"),      # Black axis lines
    axis.title = element_text(size = 14)
  )

ggsave('cdf4KS.tiff',plot = cdf_ks,dpi = 300)
###
###
wilcox.test(caical$`%G3+C3`,caical$`%G2+C2`)
t.test(caical$`%G3+C3`,caical$`%G2+C2`)

shapiro_wilk(caical$`%G1+C1`)

round(mean(caical$`%G3s+C3s`),2)
min(cui$Nc)
range(cui$CAI)
mean(cui$CAI)
###

mean(df_plot$PR2_x)
sd(df_plot$PR2_x)
