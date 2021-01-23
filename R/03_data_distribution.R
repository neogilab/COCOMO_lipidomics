##########################################################
# Script for normality check/data distribution
# This script illustrates the distribution of data
# and check normality of data
##########################################################

# Clear workspace
rm(list = ls())

# Load packages
library(tidyverse)
library(patchwork)
library(Publish) 

# Load data 
lipidomics_data <- read_csv("data/02_lipidomics_data_tidy.csv")
lipidomics_info <- read_csv("data/02_lipidomics_data_info.csv")
clinical_data <- read_csv("data/02_clinical_data_tidy.csv")



# Normality check: Density plot --------------------------------------------------------------------
# Pivot data frame 
pivot_lipidomics_data <- lipidomics_data %>% 
  pivot_longer(cols = -c(Condition, GENDER, ID_Condition), 
               names_to = "Lipid", 
               values_to = "Concentrations")

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



# Normality test: Kolmogorov-Smirnov Test ---------------------------------------------
# Calculate mean of every lipid across samples
mean_concentration <- pivot_lipidomics_data %>% 
  group_by(Lipid) %>% 
  summarize(Mean = mean(Concentrations))

# Distribution of mean lipid concentration
mean_concentration %>% 
  ggplot(aes(x = Mean)) + 
  geom_density() + 
  labs(subtitle = "Distribution of lipid concentrations",
       x = "Lipid concentrations")

# Kolmogorov-Smirnov Test check if the concentrations are normaly distributed
pval_ks_test <- ks.test(x = mean_concentration$Mean, 
        y = "pnorm", 
        mean = mean(mean_concentration$Mean), 
        sd = sd(mean_concentration$Mean))$p.value #p-value = 2.88e-13
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

# Count number of lipids in lipid categories
count_lip_categories <- lipidomics_info_extra %>% 
  group_by(`Lipid catagories`) %>% 
  summarise(n = n()) %>% 
  mutate(fraction = n/sum(n),
         percentages = round(fraction*100,2),
         ymax = cumsum(fraction),
         ymin = c(0,head(ymax, n=-1)),
         labelPosition = (ymax + ymin)/2,
         label = paste0("n = ", n,",\n", percentages, "%"))

# Donut chart of lipid categories
lipid_categories <- count_lip_categories %>% 
  ggplot(aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=`Lipid catagories`)) +
  geom_rect() +
  coord_polar(theta = "y") +
  xlim(c(2,4)) + 
  theme_void() +
  geom_label(x=3.6, aes(y=labelPosition, label=label), size = 2.5) + 
  scale_fill_brewer(palette = 2, direction = -1)
ggsave("results/03_donut_plot_lipid_categories.png", plot = lipid_categories, device = "png", width = 6, height = 6)



# Distribution of clinical variables - Compare GENDER -----------------------------------------------------------------
# Distribution of VAT stratified on condition
clinical_data %>% 
  ggplot(mapping = aes(x = VAT, fill = CONDITION)) + 
  geom_density(alpha = 0.5) + 
  labs(title = 'Distribution of "VAT" stratified on Condition', x = 'Measure of "VAT"') + 
  scale_fill_discrete(name='Condition',
                      breaks=c('HIV_NoMetS', 'HIV_MetS'),
                      labels=c('HIV without MetS', 'HIV with MetS')) + 
  facet_grid(~ GENDER)

# Distribution of CD4 stratified on condition
clinical_data %>% 
  ggplot(mapping = aes(x = VAT, fill = CONDITION)) + 
  geom_density(alpha = 0.5) + 
  labs(title = 'Distribution of "VAT" stratified on Condition', x = 'Measure of "VAT"') + 
  scale_fill_discrete(name='Condition',
                      breaks=c('HIV_NoMetS', 'HIV_MetS'),
                      labels=c('HIV without MetS', 'HIV with MetS')) + 
  facet_grid(~ GENDER)

# Distribution of CD8 stratified on condition
clinical_data %>% 
  ggplot(mapping = aes(x = SAT, fill = CONDITION)) + 
  geom_density(alpha = 0.5) + 
  labs(title = 'Distribution of "SAT" stratified on Condition', x = 'Measure of "SAT"') + 
  scale_fill_discrete(name='Condition',
                      breaks=c('HIV_NoMetS', 'HIV_MetS'),
                      labels=c('HIV without MetS', 'HIV with MetS')) +
  facet_grid(~ GENDER)

# Distribution of thymedine_exposure stratified on Condition
clinical_data %>% 
  ggplot(mapping = aes(x = as.character(thym_exposure), fill = CONDITION)) + 
  geom_histogram(stat = "count", position = position_dodge()) +
  labs(title = 'Distribution of "Thymedine exposure" stratified on "Condition" divided in "Gender"', x = 'Thymidine exposure False(0) and True(1)') + 
  scale_fill_discrete(name='Condition',
                      breaks=c('HIV_NoMetS', 'HIV_MetS'),
                      labels=c('HIV without MetS', 'HIV with MetS')) +
  facet_grid(~ GENDER)



# Distribution of clinical variables - Compare CONDITION -----------------------------------------------------------------
# Distribution of VAT stratified on thymidine exposure
clinical_data %>% 
  ggplot(mapping = aes(x = VAT, fill = as.character(thym_exposure))) + 
  geom_density(alpha = 0.5) +
  labs(title = 'Distribution of "VAT" stratified on thymidine exposure', x = 'Measure of "VAT"') +
  scale_fill_discrete(name='Thymidine Exposure',
                      breaks=c('0', '1'),
                      labels=c('False', 'True')) +
  facet_grid(~ CONDITION)

# Distribution of VAT stratified on ART NNRTI use
clinical_data %>% 
  ggplot(mapping = aes(x = VAT, fill = as.character(NNRTI))) + 
  geom_density(alpha = 0.5) +
  labs(title = 'Distribution of "VAT" stratified on ART NNRTI use', x = 'Measure of "VAT"') +
  scale_fill_discrete(name='ART NNRTI',
                      breaks=c('0','1'),
                      labels=c('No', 'Yes')) +
  facet_grid(~ CONDITION)

# Boxplot of VAT
clinical_data %>% 
  ggplot(aes(GENDER, VAT, fill = CONDITION)) +
  geom_boxplot() 

# Boxplot of AGE
clinical_data %>% 
  ggplot(aes(GENDER, AGE, fill = CONDITION)) +
  geom_boxplot() 

# Boxplot of "TAG(54:4)-FA(20:3)"
lipidomics_data %>% 
  filter(!Condition == "Ctrl") %>% 
  ggplot(aes(GENDER, `TAG(54:4)-FA(16:0)`, fill = Condition)) +
  geom_boxplot() 



# Summary statistics of clinical variables, difference between HIV_NoMetS vs. HIV_MetS -----------------------------------------------------------------
# Sex
by(clinical_data$GENDER, clinical_data$CONDITION, summary)
# Chi-squared test
chisq.test(clinical_data$GENDER, clinical_data$CONDITION)

# Age
by(clinical_data$AGE, clinical_data$CONDITION, summary)
# The data is normal if the p-value > 0.05. So the concentration variable is NOT normally distributed.
shapiro.test(clinical_data$AGE)$p.value # Not normal! 
# Left skewed distribution --> Mann Whitney
wilcox.test(AGE ~ CONDITION, data = clinical_data)

# Ethnic
by(as.character(clinical_data$Ethnic), clinical_data$CONDITION, summary)
# Chi-squared test
chisq.test(as.character(clinical_data$Ethnic), clinical_data$CONDITION)

# Immunodeficency
by(as.character(clinical_data$immunodeficiency), clinical_data$CONDITION, summary)
# Chi-squared test
chisq.test(as.character(clinical_data$immunodeficiency), clinical_data$CONDITION)

# Thymidine exposure
by(as.character(clinical_data$thym_exposure), clinical_data$CONDITION, summary)
# Chi-squared test
chisq.test(as.character(clinical_data$thym_exposure), clinical_data$CONDITION)

# VAT
by(clinical_data$VAT, clinical_data$CONDITION, summary)
# The data is normal if the p-value > 0.05. So the VAT variable is NOT normally distributed.
shapiro.test(clinical_data$VAT)$p.value # Not normal! 
# Left skewed distribution --> Mann Whitney
wilcox.test(VAT ~ CONDITION, data = clinical_data)

# SAT
by(clinical_data$SAT, clinical_data$CONDITION, summary)
# The data is normal if the p-value > 0.05. So the VAT variable is NOT normally distributed.
shapiro.test(clinical_data$SAT)$p.value # Not normal! 
# Left skewed distribution --> Mann Whitney
wilcox.test(SAT ~ CONDITION, data = clinical_data)

# NRTI
by(as.character(clinical_data$NRTI), clinical_data$CONDITION, summary)
# Chi-squared test
chisq.test(as.character(clinical_data$NRTI), clinical_data$CONDITION)

# NNRTI
by(as.character(clinical_data$NNRTI), clinical_data$CONDITION, summary)
# Chi-squared test
chisq.test(as.character(clinical_data$NNRTI), clinical_data$CONDITION)

# PI
by(as.character(clinical_data$PI), clinical_data$CONDITION, summary)
# Chi-squared test
chisq.test(as.character(clinical_data$PI), clinical_data$CONDITION)

# INSTI
by(as.character(clinical_data$INSTI), clinical_data$CONDITION, summary)
# Chi-squared test
chisq.test(as.character(clinical_data$INSTI), clinical_data$CONDITION)

# Other/Unknown ART
by(as.character(clinical_data$Other_Unknown), clinical_data$CONDITION, summary)
# Chi-squared test
chisq.test(as.character(clinical_data$Other_Unknown), clinical_data$CONDITION)



# Utables -----------------------------------------------------------------
# Change class from numeric to chracther
cols = c('Ethnic', 'immunodeficiency', 'thym_exposure','NRTI', 'NNRTI', 'PI', 'INSTI', 'Other_Unknown');    
clinical_data[,cols] = apply(clinical_data[,cols], 2, function(x) as.character(as.numeric(x)));

# Check class type
sapply(clinical_data, class)

# Use utable to characterise data
univariateTable(CONDITION ~ 
                  GENDER+AGE+Ethnic+immunodeficiency+thym_exposure+VAT+SAT+
                  NRTI+NNRTI+PI+INSTI+Other_Unknown, 
                data = clinical_data)