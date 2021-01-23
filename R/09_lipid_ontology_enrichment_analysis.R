##########################################################
# Script handles output data from the lipid ontology enrichment
# analysis made via LION/web, "target-list mode" [http://www.lipidontology.com]
# Output data consists of lipid ontology data from the 3 
# communities/modules found in the network analysis (jupyter notebook)
##########################################################

# Clear workspace
rm(list = ls())

# Load packages
library(tidyverse)
library(patchwork)

# Load data 
c1_enrich <- read_csv("data/LIONweb/c1/LION-enrichment-job1.csv")
c2_enrich <- read_csv("data/LIONweb/c2/LION-enrichment-job1.csv")
c3_enrich <- read_csv("data/LIONweb/c3/LION-enrichment-job1.csv")
network_table <- read_csv("data/08_cytoscape_pos_table_FDR.csv")

# Variables
log_thres  <- -log10(0.05)
cut_off    <- 0.7



# Bar plots --------------------------------------------------------------------
c1_plot <- c1_enrich %>%
  mutate(`Log10 adjusted p-value (FDR)` = -log10(`FDR q-value`)) %>%
  filter(`Log10 adjusted p-value (FDR)` > cut_off) %>% 
  ggplot(aes(x = reorder(Discription,`Log10 adjusted p-value (FDR)`), 
             y = `Log10 adjusted p-value (FDR)`, 
             fill = `Log10 adjusted p-value (FDR)`)) +
  geom_bar(stat = "identity", alpha = .9, width = .4) +
  scale_fill_gradient(low = "grey", high = "skyblue") +
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
  scale_fill_gradient(low = "grey", high = "red") +
  coord_flip() +
  geom_abline(slope=0, intercept = log_thres,  col = "grey", lty = 1.5) +
  xlab("") +
  ylab("Log10 adjusted p-value (FDR)") +
  ylim(0, 30) +
  theme_classic() +
  labs(fill = "Log10 adjusted\np-value (FDR)")

# Save enrichment plots
ggsave("results/09_1_LION_enrichmentplot_c1_c2.png", plot = c1_plot / c2_plot , device = "png")
ggsave("results/09_2_LION_enrichmentplot_c3.png"   , plot = c3_plot , device = "png", height = 10)



# Node degree in the network --------------------------------------------------
avg_degree <- network_table %>% 
  group_by(community) %>% 
  summarise_at(vars(-c(Biochemicals, logFC_limma, adj_pval_limma, Group, Super_pathway, sign_lip_met)), funs(round(mean(., na.rm=TRUE), 2)))

individual_degree <- network_table %>% 
  select(Biochemicals, degree, community) %>% 
  arrange(desc(degree))