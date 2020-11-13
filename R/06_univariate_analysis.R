##########################################################
# Script for univariate analysis of lipidomics data
# This script runs a non-parametric rank correlation test and calculates the fold change
# Furthermore visualizes data with a volcano plot
##########################################################

# Clear workspace
rm(list = ls())

# Load packages
library(tidyverse)
library(readr)
library(ggrepel)
library(patchwork)
library(limma)


# Variables ----------------------------------------------------
# Mann Whitney: Thresholds for the log fold change and the adjusted pvalue
log_thres <- 0.8
pval_thres <- 0.001

# Limma: Thresholds for the log fold change and the adjusted pvalue
log_thres_limma <- 0.8
pval_thres_limma <- 0.001



# Load data ---------------------------------------------------------------
lipidomics_data <- read_csv("data/02_lipidomics_data_tidy.csv")
lipidomics_info <- read_csv("data/02_lipidomics_data_info.csv")

##met_lip_data <- read_csv("data/02_metabolomics_lipidomics_data_tidy.csv")
##met_lip_data_kegg_hmdb <- read_csv("data/02_metabolomics_lipidomics_data_info.csv")



# Condition groups --------------------------------------------------------
HIV_NoMetS <- lipidomics_data %>% 
  filter(str_detect(ID_Condition, "HIV_NoMetS")) #%>% 
  #filter(GENDER == "Male") 
  
HIV_MetS <- lipidomics_data %>% 
  filter(str_detect(ID_Condition, "HIV_MetS")) #%>% 
  #filter(GENDER == "Male")

# Filter data set to contain measures for HIVnomets and HIVmets
HIVnomets_HIVmets <- lipidomics_data
  #filter(GENDER == "Male")



# Mann Whitney: HIVnomets vs. HIVmets ------------------------------------------------------------
# Check if the abundance of the lipids between HIVnomets and HIVmets differs using a Mann-Withney U test, by looking at the pvalue
pval_HIVnomets_HIVmets <- sapply(HIVnomets_HIVmets[,4:ncol(HIVnomets_HIVmets)], function(x) wilcox.test(x ~ HIVnomets_HIVmets$Condition)$p.value) %>% 
  as.data.frame() 

# Abundance ratio (log base 2 fold change)
fc_HIVnomets_HIVmets <- apply(HIV_MetS[,4:ncol(HIVnomets_HIVmets)], 2, FUN=mean)/apply(HIV_NoMetS[,4:ncol(HIVnomets_HIVmets)], 2, FUN=mean)
fclog2_HIVnomets_HIVmets <- log2(fc_HIVnomets_HIVmets)

# Correct for multiple testing by Benjamini-Hochberg to determine FDR
pval_BH_HIVnomets_HIVmets <- p.adjust(pval_HIVnomets_HIVmets$., method = "BH")

# Plot the multiple corrected pvalues, with a significance level of an FDR adjusted FDR pvalue (qvalue) < 0.001
plot(sort(pval_BH_HIVnomets_HIVmets), log = "y", pch = 16, xlab = "Lipids", ylab = "p-values", ylim = c(1e-17,1))
abline(h = 0.001, col = "red", lty = 2)

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



# Volcano plot of Mann Whitney stats: HIVnomets vs. HIVmets ------------------------------------------------------------
# Subset data matching the threshold values
#HIVnomets_HIVmets_sub = subset(stats_HIVnomets_HIVmets,abs(fclog2_HIVnomets_HIVmets) >= log_thres | pval_BH_HIVnomets_HIVmets < pval_thres)

# Volcano plot with ggplot illustrating the significant abundance between HIVnomets vs. HIVmets
mw_title <- substitute(paste("Thresholds: q-value < ", pval_thres, " and log2-value > ", log_thres, ". (HIV_NoMetS vs. HIV_MetS)", sep=""), list(pval_thres = pval_thres, log_thres = log_thres))
volcano_mw <- ggplot(stats_HIVnomets_HIVmets, aes(x = fclog2_HIVnomets_HIVmets, 
                                    y = -log10(pval_BH_HIVnomets_HIVmets),
                                    color = ifelse(abs(fclog2_HIVnomets_HIVmets)>log_thres,"grey","forestgreen"))) +
  geom_point() +
  xlim(-2, 2) +
  ylim(0, 16) + 
  xlab(expression("Fold Change, Log"[2]*"")) +  #(HIV_NoMetS vs. HIV_MetS)")) +
  ylab(expression("Benjamini-Hochberg adjusted p-value, Log"[10]*"")) +
  geom_vline(
    xintercept = c(-log_thres,log_thres),
    col = "forestgreen",
    linetype = "dotted",
    size = 1) +
  geom_hline(
    yintercept = c(-log10(pval_thres),-log10(pval_thres)),
    col = "forestgreen",
    linetype = "dotted",
    size = 1) +
  theme(legend.position = "none")+
  scale_colour_manual(values = c("grey", "forestgreen")) +
  #geom_text_repel(HIVnomets_HIVmets_sub, 
   #               mapping = aes(fclog2_HIVnomets_HIVmets, -log10(pval_BH_HIVnomets_HIVmets), label = row.names(HIVnomets_HIVmets_sub)),
    #              nudge_x      = 0.2,
     #             direction    = "y",
      #            hjust        = 0.5,
       #           segment.size = 0.1,
        #          segment.color = "grey50", 
         #         force = 1, 
          #        size = 3) +
  labs(title = "Volcano plot: Mann Whitney U test", 
       subtitle = mw_title) 

ggsave("results/06_volcano_HIV_MW.png", plot = volcano_mw, device = "png") #, width = 6.17, height = 3.1)



# Limma test (lipidomics) -------------------------------------------------------------------
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

# Fit model and compute t-statistics
fit <- lmFit(lipid_df, design)
fit <- eBayes(fit)

# Sorting by raw pvalue between lipid concentration from HIV_NoMets vs. HIV_MetS
stats_HIV_limma <- topTable(fit, sort.by = "P", n = Inf)

# Write toptable to file
stats_HIV_limma_df <- rownames_to_column(stats_HIV_limma, "Biochemicals")
stats_HIV_limma_df <- merge(x = stats_HIV_limma_df, 
                              y = as.data.frame(lipidomics_data_kegg_hmdb))
write_csv(stats_HIV_limma_df, "data/06_limmatest_HIV.csv")



# Volcano plot of limma test (lipidomics): HIVnomets vs. HIVmets ------------------------------------------------------------
# Subset data matching the threshold values
#limma_HIV_sub = subset(stats_HIV_limma,abs(stats_HIV_limma$logFC) >= log_thres_limma | stats_HIV_limma$adj.P.Val < pval_thres)

# Volcano plot with ggplot illustrating the significant abundance between HIVnomets vs. HIVmets
limma_title <- substitute(paste("Thresholds: q-value < ", pval_thres_limma, " and log2-value > ", log_thres_limma, ". (HIV_NoMetS vs. HIV_MetS)", sep=""), list(pval_thres_limma = pval_thres_limma, log_thres_limma = log_thres_limma))
volcano_limma <- ggplot(stats_HIV_limma, aes(x = logFC,
                                             y = -log10(adj.P.Val),
                                             color = ifelse(abs(logFC)>log_thres_limma,"grey","forestgreen"))) +
  geom_point() +
  xlim(-2, 2) +
  ylim(0, 16) + 
  xlab(expression("Fold Change, Log"[2]*"")) +  #(HIV_NoMetS vs. HIV_MetS)")) +
  ylab(expression("Benjamini-Hochberg adjusted p-value, Log"[10]*"")) +
  geom_vline(
    xintercept = c(-log_thres_limma,log_thres_limma),
    col = "forestgreen",
    linetype = "dotted",
    size = 1) +
  geom_hline(
    yintercept = c(-log10(pval_thres_limma),-log10(pval_thres_limma)),
    col = "forestgreen",
    linetype = "dotted",
    size = 1) +
  theme(legend.position = "none")+
  scale_colour_manual(values = c("grey", "forestgreen")) +
  #geom_text_repel(HIVnomets_HIVmets_sub, 
  #               mapping = aes(fclog2_HIVnomets_HIVmets, -log10(pval_BH_HIVnomets_HIVmets), label = row.names(HIVnomets_HIVmets_sub)),
  #              nudge_x      = 0.2,
  #             direction    = "y",
  #            hjust        = 0.5,
  #           segment.size = 0.1,
  #          segment.color = "grey50", 
  #         force = 1, 
  #        size = 3) +
  labs(title = "Volcano plot: limma test", 
       subtitle = limma_title) 
ggsave("results/06_volcano_HIV_limma.png", plot = volcano_limma, device = "png") #, width = 6.17, height = 3.1)


# Combine volcano plots Limma test + Mann whitney test (lipidomics)
ggsave("results/06_volcano_plots.png", plot = volcano_mw + volcano_limma, device = "png", width = 15) #, width = 6.17, height = 3.1)



# Significant lipids between HIV groups (lipidomics) ------------------------------------------------------
sign_lipids_MW <- subset(stats_HIV_MW, stats_HIV_MW$adj_pvalue <= pval_thres)
sign_lipids_limma <- subset(stats_HIV_limma_df, stats_HIV_limma_df$adj.P.Val <= pval_thres_limma)

method_list_univariate <- list('Mann Whitney U Test' = sign_lipids_MW$Biochemicals, 
                               'Limma test with "limma"' = sign_lipids_limma$Biochemicals)

# Save list to a file
save(method_list_univariate, file="data/06_methods_univariate.RData")









# Limma test (metabolomics + lipidomics) -------------------------------------------------------------------
#met_lip_df <- HIVnomets_HIVmets_lip_met %>% 
 # select(everything(), -c(Gender, Condition, ID_Condition)) %>% 
 # log2() %>% 
 # t()

#colnames(met_lip_df) <- HIVnomets_HIVmets_lip_met$ID_Condition
#met_lip_df <- as.data.frame(met_lip_df)

# Two groups in the lipidomics data, the control group is excluded
#group <- HIVnomets_HIVmets_lip_met %>% 
 # mutate(HIV_binary = case_when(Condition == "HIV_NoMetS" ~ "0",
    #                            Condition == "HIV_MetS" ~ "1", 
    #                            TRUE ~ Condition))
#design <- cbind(Group = as.numeric(group$HIV_binary))

# Fit model and compute t-statistics
#fit <- lmFit(met_lip_df, design)
#fit <- eBayes(fit)

# Sorting by raw pvalue between lipid concentration from HIV_NoMets vs. HIV_MetS
#stats_HIV_limma_met_lip <- topTable(fit, sort.by = "P", n = Inf)

# Write toptable to file
#stats_HIV_limma_met_lip <- rownames_to_column(stats_HIV_limma_met_lip, "Biochemicals")
#stats_HIV_limma_met_lip <- merge(x = stats_HIV_limma_met_lip, 
 #                                y = as.data.frame(met_lip_data_kegg_hmdb),
 #                                by = "Biochemicals")
#write_csv(stats_HIV_limma_met_lip, "data/06_limmatest_HIV_met_lip.csv")