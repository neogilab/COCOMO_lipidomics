##########################################################
# Script for normality check/data distribution
# This script illustrates the distribution of data
# and it runs two normality tests.
##########################################################

# Clear workspace
rm(list = ls())

# Load packages
library(tidyverse)
library(patchwork)
library(Publish) # For univariate tables

# Load data 
lipidomics_data <- read_csv("data/02_lipidomics_data_tidy.csv")
lipidomics_info <- read_csv("data/02_lipidomics_data_info.csv")

#metabolite_sign <- read_csv("data/02_metabolomics_sign_metabolites.csv")
#clinical_data <- read_csv("data/02_clinical_data_tidy.csv")
#network_table <- read_csv("data/08_cytoscape_pos_table_FDR.csv") %>% rename(Biochemicals = X1)


# Pivot data frame -------------------------------------------------
pivot_lipidomics_data <- lipidomics_data %>% 
  pivot_longer(cols = -c(Condition, GENDER, ID_Condition), 
               names_to = "Lipid", 
               values_to = "Concentrations")



# Normality check: Density plot --------------------------------------------------------------------
# Distribution on raw data values
raw_distribution <- pivot_lipidomics_data %>% 
  ggplot(aes(x = Concentrations)) + 
  geom_density() + 
  labs(subtitle = "Distribution of lipid concentrations",
       x = "Lipid concentrations")

# Distribution on log base 2 transformed data values
log2_distribution <- pivot_lipidomics_data %>% 
  ggplot(aes(x = log2(Concentrations))) + 
  geom_density() + 
  labs(subtitle = "Distribution of lipid concentrations",
       x = expression("Lipid concentrations, Log"[2]*""))

# Save in one plot
ggsave("results/03_distribution_plots.png", plot = raw_distribution + log2_distribution, device = "png", width = 6.17, height = 3.1)



# Normality test: Shapiro-Wilk Test ---------------------------------------
# Calculate mean of every lipid across samples
mean_concentration <- pivot_lipidomics_data %>% 
  group_by(Lipid) %>% 
  summarize(Mean = mean(Concentrations))

# Shapiro test, check if the concentrations are normaly distributed
pval_shapiro_test <- shapiro.test(mean_concentration$Mean)$p.value #p-value = 7.649732e-33
#The data is normal if the p-value > 0.05. So the concentration variable is NOT normally distributed.

# Distribution of mean lipid concentration
mean_concentration %>% 
  ggplot(aes(x = Mean)) + 
  geom_density() + 
  labs(subtitle = "Distribution of lipid concentrations",
       x = "Lipid concentrations")



# Normality test: Kolmogorov-Smirnov Test ---------------------------------------------
pval_ks_test <- ks.test(x = mean_concentration$Mean, 
        y = "pnorm", 
        mean = mean(mean_concentration$Mean), 
        sd = sd(mean_concentration$Mean))$p.value #p-value = 1.014642e-08
#The data is normal if the p-value > 0.05. So the concentration variable is NOT normally distributed.





# Distribution of lipid classes -------------------------------------------
lipidomics_info_extra <- lipidomics_info %>% 
  mutate(`Lipid catagories` = case_when(str_detect(Sub_pathway, "CE Ester") ~ "Sterol lipids",
                                        str_detect(Sub_pathway, "Ceramide") ~ "Sphingolipids",
                                        str_detect(Sub_pathway, "DAG Ester") ~ "Glycerolipids",
                                        str_detect(Sub_pathway, "Dihydroceramide") ~ "Sphingolipids",
                                        str_detect(Sub_pathway, "Hexosylceramide") ~ "Sphingolipids",
                                        str_detect(Sub_pathway, "Lactosylceramide") ~ "Sphingolipids",
                                        str_detect(Sub_pathway, "LPC Ester") ~ "Glycerophospholipids",
                                        str_detect(Sub_pathway, "LPE Ester") ~ "Glycerophospholipids",
                                        str_detect(Super_pathway, "MAG") ~ "Glycerolipids",
                                        str_detect(Sub_pathway, "PC Ester") ~ "Glycerophospholipids",
                                        str_detect(Sub_pathway, "PE Ester") ~ "Glycerophospholipids",
                                        str_detect(Sub_pathway, "PE Ether") ~ "Glycerophospholipids",
                                        str_detect(Sub_pathway, "PE Plasmalogen") ~ "Glycerophospholipids",
                                        str_detect(Sub_pathway, "PI Ester") ~ "Glycerophospholipids",
                                        str_detect(Sub_pathway, "Sphingomyelin") ~ "Sphingolipids",
                                        str_detect(Sub_pathway, "TAG Ester") ~ "Glycerolipids")) 


# Lipid categories
count_lip_categories <- lipidomics_info_extra %>% 
  group_by(`Lipid catagories`) %>% 
  summarise(n = n()) %>% 
  mutate(fraction = n/sum(n),
         percentages = round(fraction*100,2),
         ymax = cumsum(fraction),
         ymin = c(0,head(ymax, n=-1)),
         labelPosition = (ymax + ymin)/2,
         label = paste0("n = ", n,",\n", percentages, "%"))

# Donut chart
lipid_categories <- count_lip_categories %>% 
  ggplot(aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=`Lipid catagories`)) +
  geom_rect() +
  coord_polar(theta = "y") +
  xlim(c(2,4)) + 
  theme_void() +
  geom_label(x=3.6, aes(y=labelPosition, label=label), size = 2.5) + 
  scale_fill_brewer(palette = 2, direction = -1)
lipid_categories
ggsave("results/03_donut_plot_lipid_categories.png", plot = lipid_categories, device = "png")#, width = 6.17, height = 3.1)
