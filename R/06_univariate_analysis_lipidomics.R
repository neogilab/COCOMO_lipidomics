##########################################################
# Script for univariate analysis of lipidomics data
# This script runs a non-parametric rank correlation test
# Furthermore visualizes data with a volcano plot
##########################################################

# Clear workspace
rm(list = ls())

# Load packages
library(tidyverse)
library(patchwork)
library(limma)

# Load data 
lipidomics_data <- read_csv("data/02_lipidomics_data_tidy.csv")
lipidomics_info <- read_csv("data/02_lipidomics_data_info.csv")

# Variables: Mann Whitney & Limma thresholds (log fold change and adjusted pvalue)
log_thres_MW <- 0 
pval_thres_MW <- 0.001
log_thres_limma <- 0
pval_thres_limma <- 0.001



# Pre-filter data for HIVnomets and HIVmets -------------------------------
HIV_NoMetS <- lipidomics_data %>% 
  filter(str_detect(ID_Condition, "HIV_NoMetS")) 
HIV_MetS <- lipidomics_data %>% 
  filter(str_detect(ID_Condition, "HIV_MetS")) 
HIVnomets_HIVmets <- lipidomics_data  %>% 
  filter(!str_detect(Condition, "Ctrl"))



# Mann Whitney ------------------------------------------------------------
# Check differnce between lipid abundances between HIV-NoMetS vs. HIVMetS
pval_HIVnomets_HIVmets <- sapply(HIVnomets_HIVmets[,4:ncol(HIVnomets_HIVmets)], function(x) wilcox.test(x ~ HIVnomets_HIVmets$Condition)$p.value) %>% 
  as.data.frame() 

# Abundance ratio (log2 fold change)
fc_HIVnomets_HIVmets <- apply(HIV_MetS[,4:ncol(HIVnomets_HIVmets)], 2, FUN=mean)/apply(HIV_NoMetS[,4:ncol(HIVnomets_HIVmets)], 2, FUN=mean)
fclog2_HIVnomets_HIVmets <- log2(fc_HIVnomets_HIVmets)

# Correct for multiple testing by Benjamini-Hochberg to determine FDR
pval_BH_HIVnomets_HIVmets <- p.adjust(pval_HIVnomets_HIVmets$., method = "BH")

# Overview of significant lipids between HIVnomets and HIVmets, ordered by FDR adjusted pvalue
stats_HIVnomets_HIVmets <- data.frame(pval_HIVnomets_HIVmets, 
                                      pval_BH_HIVnomets_HIVmets, 
                                      fclog2_HIVnomets_HIVmets)[order(pval_BH_HIVnomets_HIVmets),]

# Save statistics from Mann Withney test
stats_HIV_MW <- stats_HIVnomets_HIVmets %>% 
  mutate(Biochemicals = row.names(stats_HIVnomets_HIVmets)) %>% 
  select(., Biochemicals, pvalue = ., log2FC = fclog2_HIVnomets_HIVmets, adj_pvalue = pval_BH_HIVnomets_HIVmets)
stats_HIV_MW <- merge(x = stats_HIV_MW, 
                            y = as.data.frame(lipidomics_info))

# Write MW stats to file
write_csv(stats_HIV_MW, "data/06_MWtest_HIV.csv")



# Mann Whitney: Volcano plot ------------------------------------------------------------
# Define regulated lipid abundances
stats_HIV_MW$diffexpressed <- "No"
stats_HIV_MW$diffexpressed[stats_HIV_MW$log2FC > log_thres_MW & stats_HIV_MW$adj_pvalue<pval_thres_MW] <- "Up"
stats_HIV_MW$diffexpressed[stats_HIV_MW$log2FC < log_thres_MW & stats_HIV_MW$adj_pvalue<pval_thres_MW] <- "Down"

# Volcano plot with ggplot illustrating the significant abundance between HIVnomets vs. HIVmets
volcano_mw <- ggplot(stats_HIV_MW, aes(x = log2FC, 
                                    y = -log10(adj_pvalue),
                                    col = diffexpressed)) + 
  geom_point() +
  xlim(-1.5, 1.5) +
  ylim(0, 16) + 
  xlab(expression("Fold Change, Log"[2]*"")) +  
  ylab(expression("FDR adjusted p-value, Log"[10]*"")) +
  geom_hline(
    yintercept = c(-log10(pval_thres_MW),-log10(pval_thres_MW)),
    col = "black",
    linetype = "dotted",
    size = 1) +
  theme(legend.position = "none")+
  scale_colour_manual(values = c("red","grey", "forestgreen")) +
  labs(title = "Volcano plot: Mann Whitney U test") 
ggsave("results/06_1_volcano_HIV_MW.png", plot = volcano_mw, device = "png")



# Limma test -------------------------------------------------------------------
# Check differnce between lipid abundances between HIV-NoMetS vs. HIVMetS
lipid_df <- HIVnomets_HIVmets %>% 
  select(everything(), -c(GENDER, Condition, ID_Condition)) %>% 
  log2() %>% 
  t()

colnames(lipid_df) <- HIVnomets_HIVmets$ID_Condition
lipid_df <- as.data.frame(lipid_df)

# Two groups in the lipidomics data, the control group is excluded
group <- HIVnomets_HIVmets %>% 
  mutate(HIV_binary = case_when(Condition == "HIV_NoMetS" ~ "0",
                                Condition == "HIV_MetS" ~ "1", 
                                TRUE ~ Condition))
design <- cbind(Group = as.numeric(group$HIV_binary))

# Fit model and compute moderated t-statistics
fit <- lmFit(lipid_df, design)
fit <- eBayes(fit)

# Sorting by raw pvalue between lipid concentration from HIV_NoMets vs. HIV_MetS
stats_HIV_limma <- topTable(fit, sort.by = "P", n = Inf)

# Write toptable to file
stats_HIV_limma_df <- rownames_to_column(stats_HIV_limma, "Biochemicals")
stats_HIV_limma_df <- merge(x = stats_HIV_limma_df, 
                              y = as.data.frame(lipidomics_info))
write_csv(stats_HIV_limma_df, "data/06_limmatest_HIV.csv")



# Limma test: Volcano plot -------------------------------------------------------------
# Define regulated lipid abundances
stats_HIV_limma$diffexpressed <- "No"
stats_HIV_limma$diffexpressed[stats_HIV_limma$logFC > log_thres_limma & stats_HIV_limma$adj.P.Val<pval_thres_limma] <- "Up"
stats_HIV_limma$diffexpressed[stats_HIV_limma$logFC < log_thres_limma & stats_HIV_limma$adj.P.Val<pval_thres_limma] <- "Down"

# Volcano plot with ggplot illustrating the significant abundance between HIVnomets vs. HIVmets
volcano_limma <- ggplot(stats_HIV_limma, aes(x = logFC,
                                             y = -log10(adj.P.Val),
                                             col = diffexpressed)) + 
  geom_point() +
  xlim(-1.5, 1.5) +
  ylim(0, 16) + 
  xlab(expression("Fold Change, Log"[2]*"")) +
  ylab(expression("FDR adjusted p-value, Log"[10]*"")) +
  geom_hline(
    yintercept = c(-log10(pval_thres_limma),-log10(pval_thres_limma)),
    col = "black",
    linetype = "dotted",
    size = 1) +
  theme(legend.position = "right")+
  scale_colour_manual(values = c("red","grey", "forestgreen")) +
  labs(title = "Volcano plot: Limma test", col = "Significant \ndifferential \nabundance") 
ggsave("results/06_2_volcano_HIV_limma.png", plot = volcano_limma, device = "png") #, width = 6.17, height = 3.1)

# Combine volcano plots Limma test + Mann whitney test (lipidomics)
ggsave("results/06_3_volcano_plots.png", plot = volcano_mw + volcano_limma, device = "png", width = 15) #, width = 6.17, height = 3.1)



# Significant lipids ------------------------------------------------------
sign_lipids_MW <- subset(stats_HIV_MW, stats_HIV_MW$adj_pvalue <= pval_thres_MW )
sign_lipids_limma <- subset(stats_HIV_limma_df, stats_HIV_limma_df$adj.P.Val <= pval_thres_limma)

method_list_univariate <- list('Mann Whitney U Test' = sign_lipids_MW$Biochemicals, 
                               'Limma test with "limma"' = sign_lipids_limma$Biochemicals)

# Save list to a file
save(method_list_univariate, file="data/06_methods_univariate.RData")