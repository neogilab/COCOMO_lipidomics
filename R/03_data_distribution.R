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
clinical_data <- read_csv("data/02_clinical_data_tidy.csv")

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






# Cleaning of the medical variables -> distribution of the clinical variables eventually --------------------------------------------
# Subsetting the clinical variables
clinical_data_spec_var <- clinical_data %>% 
  select(ID_Condition, AGE, GENDER, CONDITION, Ethnic, CD4, CD8, WHR,
         ART1, ART2, ART3, ART4, ART5, # Current ART use
         ART1_prev, ART2_prev, ART3_prev, ART4_prev, ART5_prev, ART6_prev, ART7_prev, ART8_prev, ART9_prev, ART10_prev) %>%  # Previous ART use
  replace_na(list(Ethnic = 4)) %>% # 4 = other/unknown
  mutate(age_groupings = findInterval(AGE, c(40, 50, 60, 70))) %>% # Make age groupings 1[40-49], 2[50-59], 3[60-69], 4[70-79]
  mutate(Gender_enc = case_when(GENDER == 'Female' ~ 0,
                                GENDER == 'Male' ~ 1)) %>% 
  mutate(Condition_enc = case_when(CONDITION == 'HIV_NoMetS' ~ 0,
                                   CONDITION == 'HIV_MetS' ~ 1)) 

# Find unique terms for current ART
Art1_curr <- as_tibble(unique(clinical_data_spec_var$ART1))
Art2_curr <- as_tibble(unique(clinical_data_spec_var$ART2))
Art3_curr <- as_tibble(unique(clinical_data_spec_var$ART3))
Art4_curr <- as_tibble(unique(clinical_data_spec_var$ART4))
Art5_curr <- as_tibble(unique(clinical_data_spec_var$ART5))

# Find unique terms for previous ART 
Art1_prev <- as_tibble(unique(clinical_data_spec_var$ART1_prev))
Art2_prev <- as_tibble(unique(clinical_data_spec_var$ART2_prev))
Art3_prev <- as_tibble(unique(clinical_data_spec_var$ART3_prev))
Art4_prev <- as_tibble(unique(clinical_data_spec_var$ART4_prev))
Art5_prev <- as_tibble(unique(clinical_data_spec_var$ART5_prev))
Art6_prev <- as_tibble(unique(clinical_data_spec_var$ART6_prev))
Art7_prev <- as_tibble(unique(clinical_data_spec_var$ART7_prev))
Art8_prev <- as_tibble(unique(clinical_data_spec_var$ART8_prev))
Art9_prev <- as_tibble(unique(clinical_data_spec_var$ART9_prev))
Art10_prev <- as_tibble(unique(clinical_data_spec_var$ART10_prev))

# Find unique terms for all ART (curr + prev)
ART_all <- c(Art1_curr$value, Art2_curr$value, Art3_curr$value, Art4_curr$value, Art5_curr$value,
             Art1_prev$value, Art2_prev$value, Art3_prev$value, Art4_prev$value, Art5_prev$value,
             Art6_prev$value, Art7_prev$value, Art8_prev$value, Art9_prev$value, Art10_prev$value)
uniq_ART_all <- as_tibble(unique(ART_all))


# Based on the unique list of ART (43 terms), each term has been grouped into the following 7 groups
NRTI            <- c('Combivir', 'Complera', 'Complera Complera', 'Descovy', 'Emtriva(emtricitabin,FTC)', 'Epivir(lamivudin,3TC)', 'Epzicom', 
                     'Kivexa', 'Lamivudin', 'Trizivir', 'Truvada', 'Viread(tenofovir,TDF)', 'Ziagen(abavacir,ABC)')
NNRTI           <- c('Edurant (rilpivirine, RPV)', 'Intelence (etravirine, ETR)', 'Invirase (saquinavir, SQV)', 'Stocrin', 
                     'Sustiva (Stocrin, efavirenz, EFV)', 'Viramune XR (nevirapine, NVP)')
PI              <- c('Amprenavir', 'Kaletra (Aluvia, lopinavir/ritonavir, LPV/r)', 'Norvir (ritonavir, RTV)', 'Norvir(ritonavir,RTV)','Prezcobix (Rezolsta, darunavir + cobicistat)', 
                     'Prezista (darunavir, DRV)', 'Reyataz (atazanavir, ATV)', 'Tybost(cobicstat)', 'Viracept (nelfinavir, NFV)')
INSTI           <- c('Isentress (raltegravir)', 'Tivicay (dolutegravir)', 'Vitekta (elvitegravir)')
comb_dif_groups <- c('Atripla', 'Genvoya', 'Stribild', 'Triumeq')
other_unknown   <- c('83', 'Fuzeon (enfuvirtide, ENF)', 'Selzentry (Celsentri, maraviroc)', 
                   'Retrovir(zidovudin,AZT)') # Only one lipo_art in "current ART cols", maybe a mistake - therefore in other group
                    # Maraviroc is an entry inhibitor, however there are few of the, thus classified "other"
lipo_ART        <- c('Crixivan (indinavir, IDV)', 'Retrovir(zidovudin,AZT)', 'Videx(didanosine,ddl)', 'Videx(didanosine,ddl) Videx(didanosine)', 'Zerit(stavudin,d4T)')


# Function to encode the above current ART groups
encode_curr_art_groups <- function(art_col) {
  art_col_group = case_when(art_col %in% NRTI ~ 1, 
                            art_col %in% NNRTI ~ 2, 
                            art_col %in% PI ~ 3, 
                            art_col %in% INSTI ~ 4, 
                            art_col %in% comb_dif_groups ~ 5, 
                            art_col %in% other_unknown ~ 6)
  return(art_col_group)
}

# Creating new columns with encoded values for the ART groups and boolean column for thymidine exposure or not. 
clinical_data_spec_var_encoded <- clinical_data_spec_var %>% 
  mutate(ART1_group = encode_curr_art_groups(ART1)) %>% 
  mutate(ART2_group = encode_curr_art_groups(ART2)) %>% 
  mutate(ART3_group = encode_curr_art_groups(ART3)) %>% 
  mutate(ART4_group = encode_curr_art_groups(ART4)) %>% 
  mutate(ART5_group = encode_curr_art_groups(ART5)) %>% 
  mutate(thym_exposure = apply(.[, c('ART1_prev', 'ART2_prev', 'ART3_prev', 'ART4_prev', 'ART5_prev', 'ART6_prev', 'ART7_prev', 'ART8_prev', 'ART9_prev', 'ART10_prev')], 
                               1, function(r) any(r %in% lipo_ART))) %>% # If a value from lipo_ART is in any of the 10 cols, then TRUE otherwise FALSE
  mutate(thym_exposure_enc = case_when(thym_exposure == 'FALSE' ~ 0,
                                thym_exposure == 'TRUE' ~ 1)) %>% 
  select(ID_Condition, CONDITION, Condition_enc, AGE, age_groupings, GENDER, Gender_enc, Ethnic, CD4, CD8, WHR, thym_exposure_enc, 
         ART1_group, ART2_group, ART3_group, ART4_group, ART5_group, ART1, ART2, ART3, ART4, ART5)




# Distribution of clinical data -------------------------------------------
######################### Compare GENDER ######################### 
# Distribution of WHR stratified on condition
clinical_data_spec_var_encoded %>% 
  ggplot(mapping = aes(x = WHR, fill = CONDITION)) + 
    geom_density(alpha = 0.5) + 
    labs(title = 'Distribution of "Waist-Hip ratio" stratified on Condition', x = 'Measure of "Waist to hip ratio"') + 
    scale_fill_discrete(name='Condition',
                        breaks=c('HIV_NoMetS', 'HIV_MetS'),
                        labels=c('HIV without MetS', 'HIV with MetS')) + 
  facet_grid(~ GENDER)

# Distribution of CD4 stratified on condition
clinical_data_spec_var_encoded %>% 
  ggplot(mapping = aes(x = CD4, fill = CONDITION)) + 
  geom_density(alpha = 0.5) + 
  labs(title = 'Distribution of "CD4 count" stratified on Condition', x = 'Measure of "Waist to hip ratio"') + 
  scale_fill_discrete(name='Condition',
                      breaks=c('HIV_NoMetS', 'HIV_MetS'),
                      labels=c('HIV without MetS', 'HIV with MetS')) + 
  facet_grid(~ GENDER)

# Distribution of CD8 stratified on condition
clinical_data_spec_var_encoded %>% 
  ggplot(mapping = aes(x = CD8, fill = CONDITION)) + 
  geom_density(alpha = 0.5) + 
  labs(title = 'Distribution of "CD8 count" stratified on Condition', x = 'Measure of "Waist to hip ratio"') + 
  scale_fill_discrete(name='Condition',
                      breaks=c('HIV_NoMetS', 'HIV_MetS'),
                      labels=c('HIV without MetS', 'HIV with MetS')) +
  facet_grid(~ GENDER)

# Distribution of thymedine_exposure stratified on Condition
clinical_data_spec_var_encoded %>% 
  ggplot(mapping = aes(x = as.character(thym_exposure_enc), fill = CONDITION)) + 
  geom_histogram(stat = "count", position = position_dodge()) +
  labs(title = 'Distribution of "Thymedine exposure" stratified on "Condition" divided in "Gender"', x = 'Thymidine exposure False(0) and True(1)') + 
  scale_fill_discrete(name='Condition',
                      breaks=c('HIV_NoMetS', 'HIV_MetS'),
                      labels=c('HIV without MetS', 'HIV with MetS')) +
  facet_grid(~ GENDER)

######################### Compare CONDITION ######################### 
# Distribution of WHR stratified on thymidine exposure
clinical_data_spec_var_encoded %>% 
  ggplot(mapping = aes(x = WHR, fill = as.character(thym_exposure_enc))) + 
  geom_density(alpha = 0.5) +
  labs(title = 'Distribution of "Waist-Hip ratio" stratified on thymidine exposure', x = 'Measure of "Waist to hip ratio"') +
  scale_fill_discrete(name='Thymidine Exposure',
                      breaks=c('0', '1'),
                      labels=c('False', 'True')) +
  facet_grid(~ CONDITION)

# Distribution of WHR stratified on age ranges
clinical_data_spec_var_encoded %>% 
  ggplot(mapping = aes(x = WHR, fill = as.character(age_groupings))) + 
  geom_density(alpha = 0.5) +
  labs(title = 'Distribution of "Waist-Hip ratio" stratified on age ranges', x = 'Measure of "Waist to hip ratio"') +
  scale_fill_discrete(name='Age groups',
                      breaks=c('1', '2', '3', '4'),
                      labels=c('40-49', '50-59', '60-69', '70-79')) +
  facet_grid(~ CONDITION)

# Distribution of WHR stratified on ART1 use
clinical_data_spec_var_encoded %>% 
  ggplot(mapping = aes(x = WHR, fill = as.character(ART1_group))) + 
  geom_density(alpha = 0.5) +
  labs(title = 'Distribution of "Waist-Hip ratio" stratified on ART1 use', x = 'Measure of "Waist to hip ratio"') +
  scale_fill_discrete(name='ART groups',
                      breaks=c('1', '2', '3', '4', '5','6'),
                      labels=c('NRTI', 'NNRTI', 'PI', 'INSTI', 'cART (between diff. groups)', 'Other/unknown')) +
  facet_grid(~ CONDITION)

# Distribution of Gender stratified on age groups --> this plot could imply exclusion of females
clinical_data_spec_var_encoded %>% 
  ggplot(mapping = aes(x = GENDER, fill = as.character(age_groupings))) + 
  geom_histogram(stat = "count") +
  labs(title = 'Distribution of "Gender" stratified on "Age groups" divided in "Conditions"', x = 'Gender') +
  scale_fill_discrete(name='Age groups',
                      breaks=c('1', '2', '3', '4'),
                      labels=c('40-49', '50-59', '60-69', '70-79')) +
  facet_grid(~ CONDITION)


# Distribution of WHR stratified on ART1 use
clinical_data_spec_var_encoded %>% 
  ggplot(mapping = aes(x = WHR, fill = as.character(ART1_group))) + 
  geom_histogram(alpha = 0.5) +
  labs(title = 'Distribution of "Waist-Hip ratio" stratified on ART1 use', x = 'Measure of "Waist to hip ratio"') 



######################################################################################
######################################################################################
# Summary statistics of clinical variables, difference between HIV_NoMetS vs. HIV_MetS
######################################################################################
######################################################################################
# Age -----------------------------
# Density plot, identifying the distribution
clinical_data_spec_var_encoded %>% 
  ggplot(mapping = aes(x = AGE)) + 
  geom_density(alpha = 0.5) 

# Left skewed distribution --> Mann Whitney
wilcox.test(AGE ~ CONDITION, data = clinical_data_spec_var_encoded)

# Summary statistics
by(clinical_data_spec_var_encoded$AGE, clinical_data_spec_var_encoded$CONDITION, summary)



# Age gropings ------------------------------------------------------------------
# Summary statistics
by(as.character(clinical_data_spec_var_encoded$age_groupings), clinical_data_spec_var_encoded$CONDITION, summary)

# Chi-squared test
chisq.test(clinical_data_spec_var_encoded$age_groupings, clinical_data_spec_var_encoded$CONDITION)



# Gender ------------------------------------------------------------------
# Summary statistics
by(clinical_data_spec_var_encoded$GENDER, clinical_data_spec_var_encoded$CONDITION, summary)

# Chi-squared test
chisq.test(clinical_data_spec_var_encoded$GENDER, clinical_data_spec_var_encoded$CONDITION)



# Ethnic ------------------------------------------------------------------
# Summary statistics
by(as.character(clinical_data_spec_var_encoded$Ethnic), clinical_data_spec_var_encoded$CONDITION, summary)

# Chi-squared test
chisq.test(as.character(clinical_data_spec_var_encoded$Ethnic), clinical_data_spec_var_encoded$CONDITION)



# CD4 -----------------------------
# Density plot, identifying the distribution
clinical_data_spec_var_encoded %>% 
  ggplot(mapping = aes(x = CD4)) + 
    geom_density(alpha = 0.5) 

# Normal distribution --> T-test
t.test(CD4 ~ CONDITION, data = clinical_data_spec_var_encoded)

# Summary statistics
by(clinical_data_spec_var_encoded$CD4, clinical_data_spec_var_encoded$CONDITION, summary)



# CD8 -----------------------------
clinical_data_spec_var_encoded %>% 
  ggplot(mapping = aes(x = log2(CD8))) + 
  geom_density(alpha = 0.5)

# Left skewed distribution --> Mann Whitney
wilcox.test(CD8 ~ CONDITION, data = clinical_data_spec_var_encoded)

# Summary statistics
by(clinical_data_spec_var_encoded$CD8, clinical_data_spec_var_encoded$CONDITION, summary)



# Waist-hip ratio (WHR) -----------------------------
clinical_data_spec_var_encoded %>% 
  ggplot(mapping = aes(x = WHR)) + 
  geom_density(alpha = 0.5)

# Normal distribution --> T-test
t.test(WHR ~ CONDITION, data = clinical_data_spec_var_encoded)

# Summary statistics
by(clinical_data_spec_var_encoded$WHR, clinical_data_spec_var_encoded$CONDITION, summary)



# Thymidine exposure ------------------------------------------------------------------
# Summary statistics
by(as.character(clinical_data_spec_var_encoded$thym_exposure_enc), clinical_data_spec_var_encoded$CONDITION, summary)

# Chi-squared test
chisq.test(as.character(clinical_data_spec_var_encoded$thym_exposure_enc), clinical_data_spec_var_encoded$CONDITION)



# ART1_group ------------------------------------------------------------------
# Summary statistics
by(as.character(clinical_data_spec_var_encoded$ART1_group), clinical_data_spec_var_encoded$CONDITION, summary)

# Chi-squared test
chisq.test(as.character(clinical_data_spec_var_encoded$ART1_group), clinical_data_spec_var_encoded$CONDITION)



# ART2_group ------------------------------------------------------------------
# Summary statistics
by(as.character(clinical_data_spec_var_encoded$ART2_group), clinical_data_spec_var_encoded$CONDITION, summary)

# Chi-squared test
chisq.test(as.character(clinical_data_spec_var_encoded$ART2_group), clinical_data_spec_var_encoded$CONDITION)



# ART3_group ------------------------------------------------------------------
# Summary statistics
by(as.character(clinical_data_spec_var_encoded$ART3_group), clinical_data_spec_var_encoded$CONDITION, summary)

# Chi-squared test
chisq.test(as.character(clinical_data_spec_var_encoded$ART3_group), clinical_data_spec_var_encoded$CONDITION)



# ART4_group ------------------------------------------------------------------
# Summary statistics
by(as.character(clinical_data_spec_var_encoded$ART4_group), clinical_data_spec_var_encoded$CONDITION, summary)

# Chi-squared test
chisq.test(as.character(clinical_data_spec_var_encoded$ART4_group), clinical_data_spec_var_encoded$CONDITION)



# ART5_group ------------------------------------------------------------------
# Summary statistics
by(as.character(clinical_data_spec_var_encoded$ART5_group), clinical_data_spec_var_encoded$CONDITION, summary)

# Chi-squared test
chisq.test(as.character(clinical_data_spec_var_encoded$ART5_group), clinical_data_spec_var_encoded$CONDITION)






























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
