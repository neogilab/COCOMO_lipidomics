##########################################################
# Script handles output data from the lipid ontology enrichment
# analysis made via LION/web [ref. LION/web article]. 
# Output data consists of lipid ontology data from the 3 
# communities/modules found in the network analysis (jupyter notebook)
##########################################################

# Clear workspace
rm(list = ls())

# Load packages
library(tidyverse)
library(patchwork)


# Variables
log_thres  <- -log10(0.05)
cut_off    <- 0.7


# Load data ---------------------------------------------------------------
c1_enrich <- read_csv("data/LIONweb/c1_LION/c1_LION-enrichment.csv")
c2_enrich <- read_csv("data/LIONweb/c2_LION/c2_LION-enrichment.csv")
c3_enrich <- read_csv("data/LIONweb/c3_LION/c3_LION-enrichment.csv")


# Bar plots --------------------------------------------------------------------
c1_plot <- c1_enrich %>%
  mutate(`Log10 adjusted p-value (FDR)` = -log10(`FDR q-value`)) %>%
  filter(`Log10 adjusted p-value (FDR)` > cut_off) %>% 
  ggplot(aes(x = reorder(Discription,`Log10 adjusted p-value (FDR)`), 
             y = `Log10 adjusted p-value (FDR)`, 
             fill = `Log10 adjusted p-value (FDR)`)) +
  geom_bar(stat = "identity", alpha = .9, width = .4) +
  scale_fill_gradient(low = "grey", high = "green") +
  coord_flip() +
  geom_abline(slope=0, intercept = log_thres,  col = "grey", lty = 1.5) +
  xlab("") +
  ylab("") +
  ylim(0, 30) +
  theme_classic() +
  labs(fill = "Log10 adjusted\np-value (FDR)")


c2_plot <- c2_enrich %>%
  mutate(`Log10 adjusted p-value (FDR)` = -log10(`FDR q-value`)) %>%
  filter(`Log10 adjusted p-value (FDR)` > cut_off) %>% 
  ggplot(aes(x = reorder(Discription,`Log10 adjusted p-value (FDR)`), 
             y = `Log10 adjusted p-value (FDR)`, 
             fill = `Log10 adjusted p-value (FDR)`)) +
  geom_bar(stat = "identity", alpha = .9, width = .4) +
  scale_fill_gradient(low = "grey", high = "green") +
  coord_flip() +
  geom_abline(slope=0, intercept = log_thres,  col = "grey", lty = 1.5) +
  xlab("") +
  ylab("Log10 adjusted p-value (FDR)") +
  ylim(0, 30) + 
  theme_classic() +
  labs(fill = "Log10 adjusted\np-value (FDR)")


c3_plot <- c3_enrich %>%
  mutate(`Log10 adjusted p-value (FDR)` = -log10(`FDR q-value`)) %>%
  filter(`Log10 adjusted p-value (FDR)` > cut_off) %>% 
  ggplot(aes(x = reorder(Discription,`Log10 adjusted p-value (FDR)`), 
             y = `Log10 adjusted p-value (FDR)`, 
             fill = `Log10 adjusted p-value (FDR)`)) +
  geom_bar(stat = "identity", alpha = .9, width = .4) +
  scale_fill_gradient(low = "grey", high = "green") +
  coord_flip() +
  geom_abline(slope=0, intercept = log_thres,  col = "grey", lty = 1.5) +
  xlab("") +
  ylab("Log10 adjusted p-value (FDR)") +
  ylim(0, 30) +
  theme_classic() +
  labs(fill = "Log10 adjusted\np-value (FDR)")

  
# Save enrichment plots
ggsave("results/09_LION_enrichmentplot_c1_c2.png", plot = c1_plot / c2_plot , device = "png")
ggsave("results/09_LION_enrichmentplot_c3.png"   , plot = c3_plot , device = "png", height = 10)