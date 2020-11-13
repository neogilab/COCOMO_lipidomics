##########################################################
# Script connecting all results of significant lipids
# from the different analysis methods
##########################################################

# Clear workspace
rm(list = ls())

# Load packages
library(UpSetR)
library(RVenn)
library(tidyverse)
library(broom)

# Load lists from files (03_univariate_analysis.R and 05_supervised_analysis.R)
load("data/05_methods_supervised.RData")
load("data/06_methods_univariate.RData")

# Load lipidomics data 
lipidomics_data_tidy <- read_csv("data/02_lipidomics_data_tidy.csv")

# Concatenate imported lists
method_list <- c(method_list_univariate, method_list_supervised)

# UpsetR plot showing the intersection of lipids across the different methods
png("results/07_comparison_of_analyses.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # font size)      
upset(data = fromList(method_list), 
      nsets = 5, 
      order.by = "freq",
      matrix.color = "forestgreen", 
      main.bar.color = "forestgreen",
      sets.bar.color = "forestgreen")
dev.off()


# Consensus of significant lipids ------------------------------------------------------
method_list <- Venn(method_list)

# Venn diagram of the supervised methods
ggvenn(method_list, slice = c(1,2,3))

## Intersection of: '`Random Forest with "MUVR"`', '`PLS-DA with "ropls"`', '`Limma test with "limma"`' and '`Mann Whitney U Test`'
intersection <- overlap(method_list, c(1, 2, 3, 4))

# Save list to a file
intersection_list_sig_lipids <- list('Significant lipids' = intersection)
save(intersection_list_sig_lipids, file="data/07_significant_lipids.RData")

# Distribution boxplot of significant lipids ------------------------------------------------------------
# Pivot data frame 
pivot_lipidomics_data_tidy <- lipidomics_data_tidy %>% 
  pivot_longer(cols = -c(Condition, GENDER, ID_Condition), 
               names_to = "Lipid", 
               values_to = "Concentrations")

# Reorder the conditions order
pivot_lipidomics_data_tidy$Condition <- factor(pivot_lipidomics_data_tidy$Condition , levels=c("Ctrl", "HIV_NoMetS", "HIV_MetS"))

# Boxplot of significant lipids
sign_lipids_plot <- pivot_lipidomics_data_tidy %>% 
  filter(Lipid %in% intersection) %>% #, !Condition == "Ctrl") %>% 
  ggplot(mapping = aes(x = log2(Concentrations), fill = Condition)) + 
  geom_boxplot(alpha = 0.5) + 
  coord_flip() + 
  labs(title = 'Boxplot of lipid concentration', 
       subtitle = 'Significant lipids',
       x = 'Log2 of lipid concentrations') + 
  scale_fill_discrete(name="Condition") +   
  facet_wrap(~ Lipid) +
  scale_fill_manual(values = c("#32CD32", "#6495ED"))  # (Green, Blue)


sign_lipids_plot
ggsave("results/07_boxplot_significant_lipids.png", plot = sign_lipids_plot, device = "png")


# PCA of significant lipids ---------------------------
# Remove the Ctrl from the data set and extract the consensus lipids found above
lipidomics_data_tidy <- lipidomics_data_tidy %>% 
#  filter(!str_detect(ID_Condition, "Ctrl")) %>% 
  select(ID_Condition, Condition, GENDER, intersection) 

# Log2 transform data
lipidomics_data_tidy[, 4:ncol(lipidomics_data_tidy)] <- sapply(log2(lipidomics_data_tidy[, 4:ncol(lipidomics_data_tidy)]), as.numeric)
lipidomics_data_tidy_df <- as.data.frame(lipidomics_data_tidy)
rownames(lipidomics_data_tidy_df) <- lipidomics_data_tidy_df$ID_Condition

# Create PCA object 
data_pca <- lipidomics_data_tidy %>%
  select(everything(), -c(GENDER, Condition, ID_Condition)) %>% 
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


# Reorder the conditions order
data_pca_aug$Condition <- factor(data_pca_aug$Condition , levels=c("Ctrl", "HIV_NoMetS", "HIV_MetS"))

# Plot PCA with medical condition as labels
pca_condition <- data_pca_aug %>% 
  ggplot(aes(x = .fittedPC1,
             y = .fittedPC2,
             colour = Condition)) + 
  geom_point(size = 3, alpha = 0.5) + 
  labs(x = x, y = y, title = "PCA", 
       subtitle = "Significant lipids",
       color = "Condition") +
  scale_color_manual(values = c("#32CD32", "#6495ED")) + # (Green, Blue, Red - "#FA8072")
  stat_ellipse(level = 0.9)
pca_condition
ggsave("results/07_pca_significant_lipids.png", plot = pca_condition, device = "png", width = 6.17, height = 3.1)


# Heatmap & pairwise alignment of significant lipids ----------------------------------------------
# Look in Jupyter notebook in Python folder "07_consensus_analysis"


# Save file with significant lipids ---------------------------------------
as_tibble(intersection_list_sig_lipids) %>% 
  write_csv("data/07_sign_lipids.csv")
