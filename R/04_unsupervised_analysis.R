##########################################################
# Script for unsupervised analysis of lipidomics data
# This script runs dimensionality reduction by PCA
# Furthermore, cluster analysis is performed by K-means and
# hierarchichal clustering by a heatmap and a cluster dendrogram.
##########################################################

# Clear workspace
rm(list = ls())

# Load packages
library(tidyverse)
library(broom)
library(patchwork)
library(ropls)
library(ggrepel)
library(gplots)
library(RColorBrewer)
library(matrixStats)
library(dendextend)



# Set seed ----------------------------------------------------------------
# Set a seed in order to be able to create reproducible data. Sample seed by the function: sample(1e6, 1)
set.seed(218767)



# Load data ---------------------------------------------------------------
lipidomics_data_tidy <- read_csv("data/02_lipidomics_data_tidy.csv")  

##lipidomics_data_tidy <- lipidomics_data_tidy %>% 
  ##filter(!str_detect(GENDER, "Female")) 



# PCA ---------------------------------------------------------------------
# Create PCA object and log transform the variables
data_pca <- lipidomics_data_tidy %>%
  select(everything(), -c(GENDER, Condition, ID_Condition)) %>% 
  log2() %>% 
  prcomp(scale. = TRUE)

# Tidy data in order to get proper data table format
data_pca_tidy <- data_pca %>% 
  tidy("pcs")

# Scree plot using broom to tidy
scree <- data_pca_tidy %>% 
  filter(!percent < 0.01) %>% 
  mutate(variance = percent*100) %>% 
  ggplot(aes(x = PC, y = variance)) +
  geom_col() + 
  labs(title = "Scree", 
       x = "PC (cut-off < 1%)",
       y = "Proportion of variance (%)") + 
  ylim(0,100)

# Augment data in order to get a complete table with original values and PC values. 
data_pca_aug <- data_pca %>% 
  augment(lipidomics_data_tidy)

# Adding percentage to the PCA plot
x <- data_pca_tidy %>% 
  filter(PC == 1) %>% 
  pull(percent)
x <- str_c("PC1 (", round(x*100, 2), "%)")

y <- data_pca_tidy %>% 
  filter(PC == 2) %>% 
  pull(percent)
y <- str_c("PC2 (", round(y*100, 2), "%)")

# Subset outliers
pca_outliers_sub = subset(data_pca_aug, abs(.fittedPC1) >= 45 | abs(.fittedPC2) >= 23)

# Plot PCA with medical condition as labels
pca_condition <- data_pca_aug %>% 
  ggplot(aes(x = .fittedPC1,
             y = .fittedPC2,
             colour = Condition)) + 
  geom_point(size = 3, alpha = 0.5) + 
  labs(x = x, y = y, title = "PCA", color = "Condition") +
  stat_ellipse(level = 0.2) + 
  stat_ellipse(level = 0.6) +
  geom_text_repel(pca_outliers_sub, mapping = aes(x = .fittedPC1,
                                y = .fittedPC2,
                                label = ID_Condition),
                  nudge_x      = 0.2,
                  direction    = "y",
                  hjust        = 0.5,
                  segment.size = 0.01,
                  segment.color = "grey50", 
                  force = 1, 
                  size = 2) 
ggsave("results/04_pca_condition.png", plot = pca_condition, device = "png", width = 6.17, height = 3.1)

# Save PCA and scree plot
ggsave("results/04_pca_and_scree.png", plot = pca_condition + scree, device = "png", width = 8, height = 3.1)



# PCA (ropls package) ---------------------------
lipidomics_data_tidy[, 4:ncol(lipidomics_data_tidy)] <- sapply(log2(lipidomics_data_tidy[, 4:ncol(lipidomics_data_tidy)]), as.numeric)
lipidomics_data_tidy_df <- as.data.frame(lipidomics_data_tidy)
rownames(lipidomics_data_tidy_df) <- lipidomics_data_tidy_df$ID_Condition

# PCA with scores and outliers
png("results/04_pca_ropls.png")
lipidomics_pca <- opls(lipidomics_data_tidy_df[, 4:ncol(lipidomics_data_tidy_df)])
dev.off()

# PCA divided on condition group
png("results/04_pca_ropls_condition.png")
plot(lipidomics_pca, typeVc = "x-score",
     parAsColFcVn = lipidomics_data_tidy_df$Condition)
dev.off()



# K-means clustering -----------------------------------------------------------------
# Perform kmeans
data_kmeans <- data_pca_aug %>%
  select(contains("PC")) %>% 
  kmeans(centers = 2, nstart = 24)

# Add cluster column to augmented pca table
data_kmeans_aug <- data_kmeans %>%
  augment(data_pca_aug) %>%
  rename(Cluster = .cluster)

# Plot kmeans on two first principal components
kmeans_condition <- data_kmeans_aug %>% 
  ggplot(aes(x = .fittedPC1,
             y = .fittedPC2,
             colour = Cluster)) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = x, y = y, title = "K-means") +
  stat_ellipse(level = 0.2) + 
  stat_ellipse(level = 0.6) + 
  geom_text_repel(data_kmeans_aug, mapping = aes(x = .fittedPC1,
                                   y = .fittedPC2,
                                   label = ID_Condition),
                  nudge_x      = 0.2,
                  direction    = "y",
                  hjust        = 0.5,
                  segment.size = 0.01,
                  segment.color = "grey50", 
                  force = 1, 
                  size = 2) 
kmeans_condition
ggsave("results/04_kmeans_condition.png", plot = kmeans_condition, device = "png", width = 6.17, height = 3.1)



# PCA and K-means plot --------------------------------
# Save PCA and K-means as multi-panel plot
ggsave("results/04_kmeans_vs_pca.png", plot = pca_condition + kmeans_condition, device = "png", width = 8, height = 3.1)



# Hierachical clustering --------------------------------------------------
# Function to scale data from 0 to 1
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

# Filter control out of data
lip_data_log2 <- lipidomics_data_tidy %>% 
  filter(!Condition == "Ctrl") %>% 
  arrange(Condition)

# Scale log transformed data to range from 0 to 1
lip_data_norm <- as.matrix(lapply(lip_data_log2[,4:ncol(lip_data_log2)], min_max_norm))
lip_data_norm <- as.matrix.data.frame(lip_data_norm)

# Add rownames and colnames to data frame
rownames(lip_data_norm) <- colnames(lip_data_log2[4:ncol(lip_data_log2)])
colnames(lip_data_norm) <- lip_data_log2$Condition

# Calculate standard deviation estimate and order the sd according to decreasing value
sd <- rowSds(lip_data_norm, na.rm = TRUE)
o_row <- order(sd, decreasing = TRUE)[1:50]

# Order columns according to conditions
o_col <- order(colnames(lip_data_norm))

# Draw heatmap, which shows the rows with most standard deviation
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
png("results/04_heatmap.png",            
    width = 10*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # font size)      
data_hm <- heatmap.2(lip_data_norm[o_row,o_col] ,
          col = rev(hmcol), 
          main = "Heatmap of significant lipids", # heat map title
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          dendrogram="none",     # only draw a row dendrogram
          Colv=FALSE)
dev.off()



# Dendrogram --------------------------------------------------------------
data_hc <- hclust(dist(t(lip_data_norm)))
png("results/04_hierarchical_cluster_1.png",            
    width = 10*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 4)        # font size)  
plot(data_hc)
dev.off()


# Color code according to the conditions
groupCodes <- c(rep("HIV_MetS",100), rep("HIV_NoMetS",100))
colnames(lip_data_norm) <- make.unique(groupCodes)
colorCodes <- c(HIV_MetS="red", HIV_NoMetS="green")

# Make distance matrix and the hierarchical clustering
data_dist <- dist(t(lip_data_norm))
data_hc <- hclust(data_dist)
dend <- as.dendrogram(data_hc)

# Assigning the labels of dendrogram object with colors:
labels_colors(dend) <- colorCodes[groupCodes][order.dendrogram(dend)]

# Plotting dendrogram
png("results/04_hierarchical_cluster_2.png",            
    width = 10*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 4)        # font size) 
plot(dend)
dev.off()