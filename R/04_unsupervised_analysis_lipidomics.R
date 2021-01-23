##########################################################
# Script for unsupervised analysis of lipidomics data
# This script runs dimensionality reduction by PCA
# And cluster analysis is performed by K-means.
##########################################################

# Clear workspace
rm(list = ls())

# Load packages
library(tidyverse)
library(broom)
library(patchwork)
library(ropls)

# Load data 
lipidomics_data_tidy <- read_csv("data/02_lipidomics_data_tidy.csv")  

# Set a seed for reproducible data. Seed sampled with -> sample(1e6, 1)
set.seed(218767)



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

data_pca_aug$Condition <- factor(data_pca_aug$Condition , levels=c("Ctrl", "HIV_NoMetS", "HIV_MetS"))

# Plot PCA with medical condition as labels
pca_condition <- data_pca_aug %>% 
  ggplot(aes(x = .fittedPC1,
             y = .fittedPC2,
             colour = Condition)) + 
  geom_point(size = 3, alpha = 0.5) + 
  labs(x = x, y = y, title = "PCA", color = "Condition") +
  scale_color_manual(values = c("#32CD32", "#6495ED")) + 
  stat_ellipse(level = 0.95) # illustrating the 95% CI interval
ggsave("results/04_1_pca_condition.png", plot = pca_condition, device = "png", width = 6.17, height = 3.1)

# Save combined plot with PCA and scree plot
ggsave("results/04_2_pca_and_scree.png", plot = pca_condition + scree, device = "png", width = 8, height = 3.1)



# PCA (ropls package) ---------------------------
lipidomics_data_tidy[, 4:ncol(lipidomics_data_tidy)] <- sapply(log2(lipidomics_data_tidy[, 4:ncol(lipidomics_data_tidy)]), as.numeric)
lipidomics_data_tidy_df <- as.data.frame(lipidomics_data_tidy)
rownames(lipidomics_data_tidy_df) <- lipidomics_data_tidy_df$ID_Condition

# PCA with scores and outliers
png("results/04_3_pca_ropls.png")
lipidomics_pca <- opls(lipidomics_data_tidy_df[, 4:ncol(lipidomics_data_tidy_df)])
dev.off()

# PCA divided on Condition[HIV-NoMetS/HIV-MetS]
png("results/04_4_pca_ropls_condition.png")
plot(lipidomics_pca, typeVc = "x-score",
     parAsColFcVn = lipidomics_data_tidy_df$Condition)
dev.off()

# PCA divided on Sex[Male/Female]
png("results/04_5_pca_ropls_sex.png")
plot(lipidomics_pca, typeVc = "x-score",
     parAsColFcVn = lipidomics_data_tidy_df$GENDER)
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
  stat_ellipse(level = 0.95) 
ggsave("results/04_6_kmeans_condition.png", plot = kmeans_condition, device = "png", width = 6.17, height = 3.1)

# Save PCA and K-means as multi-panel plot
ggsave("results/04_7_kmeans_vs_pca.png", plot = pca_condition + kmeans_condition, device = "png", width = 8, height = 3.1)
