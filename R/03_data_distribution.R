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

# Load data 
lipidomics_data <- read_csv("data/02_lipidomics_data_tidy.csv")

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






# Distribution box plot log base 2 transformed ------------------------------------------------------------
# Significant lipids, move to another script!!
##sign_lipids_plot <- pivot_lipidomics_data %>% 
  ##filter(Lipid == c("TAG54:4-FA16:0"), !Condition == "Ctrl") %>% 
  ##ggplot(mapping = aes(x = log2(Concentrations), fill = Condition)) + 
  ##geom_boxplot(alpha = 0.5) + 
  ##coord_flip() + 
  ##labs(title = 'Boxplot of "TAG54:4-FA16:0" lipid concentration', 
    ##   x = 'Log2 of lipid concentrations') + 
  ##scale_fill_discrete(name="Condition") #+   
 # facet_wrap(~ Lipid) 
##ggsave("results/03_boxplot_significant_lipids.png", plot = sign_lipids_plot, device = "png", width = 6.17, height = 3.1)




# Distribution plots for lipid CE(12:00) -------------------------------------------------------------
# Distribution density plot
ggplot(data = lipidomics_data, mapping = aes(x = lipidomics_data$`CE(12:0)`, fill = lipidomics_data$Condition)) + 
  geom_density(alpha = 0.5) + 
  labs(title = 'Distribution of "CE(12:0)" between three condition groups', x = 'Abundance of "CE(12:0)"') + 
  scale_fill_discrete(name="Experimental\nCondition",
                      breaks=c("Ctrl", "HIV_NoMetS", "HIV_MetS"),
                      labels=c("Control", "HIV without MetS", "HIV with MetS"))

# Distribution density plot log base 2 transformed 
lipidomics_data %>% 
  ggplot(mapping = aes(x = log2(lipidomics_data$`CE(12:0)`), fill = lipidomics_data$Condition)) + 
  geom_density(alpha = 0.5) + 
  labs(title = 'Distribution of "CE(12:0)" between three condition groups', x = 'Log base 2 abundance of "CE(12:0)"') + 
  scale_fill_discrete(name="Experimental\nCondition",
                      breaks=c("Ctrl", "HIV_NoMetS", "HIV_MetS"),
                      labels=c("Control", "HIV without MetS", "HIV with MetS"))
  
# Distribution box plot log base 2 transformed 
lipidomics_data %>% 
  ggplot(mapping = aes(x = log2(lipidomics_data$`CE(12:0)`), fill = lipidomics_data$Condition)) + 
  geom_boxplot(alpha = 0.5) + 
  coord_trans() + 
  labs(title = 'Distribution of "CE(12:0)" between three condition groups', x = 'Log base 2 abundance of "CE(12:0)"') + 
  scale_fill_discrete(name="Experimental\nCondition",
                      breaks=c("Ctrl", "HIV_NoMetS", "HIV_MetS"),
                      labels=c("Control", "HIV without MetS", "HIV with MetS"))

# Distribution histogram log base 2 transformed 
ggplot(data = lipidomics_data, mapping = aes(x = log2(lipidomics_data$`CE(12:0)`), fill = lipidomics_data$Condition)) + 
  geom_histogram(alpha = 0.5) + 
  facet_wrap(~ Condition) + 
  labs(title = 'Distribution of "CE(12:0)" between three condition groups', x = 'Log base 2 abundance of "CE(12:0)"') + 
  scale_fill_discrete(name="Experimental\nCondition",
                      breaks=c("Ctrl", "HIV_NoMetS", "HIV_MetS"),
                      labels=c("Control", "HIV without MetS", "HIV with MetS"))