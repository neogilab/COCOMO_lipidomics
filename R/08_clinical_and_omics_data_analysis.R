##########################################################
# Script for regression analysis
##########################################################

# Clear workspace
rm(list = ls())

# Load packages
library(tidyverse)
library(patchwork)
library(Publish) # For univariate tables

# Load data 
lipidomics_data <- read_csv("data/02_lipidomics_data_tidy.csv")
metabolite_sign <- read_csv("data/02_metabolomics_sign_metabolites.csv")
clinical_data <- read_csv("data/02_clinical_data_tidy.csv")
network_table <- read_csv("data/08_cytoscape_pos_table_FDR.csv") 



# Cleaning of the medical variables -> distribution of the clinical variables eventually --------------------------------------------
# Subsetting the clinical variables
clinical_data_spec_var <- clinical_data %>% 
  select(ID_Condition, AGE, GENDER, CONDITION, Ethnic, CD4_nadir, CDCAIDS, VAT, SAT,
         ART1, ART2, ART3, ART4, ART5, 
         ART1_prev, ART2_prev, ART3_prev, ART4_prev, ART5_prev, ART6_prev, ART7_prev, ART8_prev, ART9_prev, ART10_prev) %>% 
  replace_na(list(Ethnic = 4)) %>% # 4 = other/unknown
  replace_na(list(VAT = mean(.$VAT, na.rm = TRUE))) %>% 
  replace_na(list(SAT = mean(.$SAT, na.rm = TRUE))) %>% 
  #mutate(age_groupings = findInterval(AGE, c(40, 50, 60, 70))) %>% # Make age groupings 1[40-49], 2[50-59], 3[60-69], 4[70-79]
  mutate(Gender_enc = case_when(GENDER == 'Female' ~ 0,
                                GENDER == 'Male' ~ 1)) %>% 
  mutate(Condition_enc = case_when(CONDITION == 'HIV_NoMetS' ~ 0,
                                   CONDITION == 'HIV_MetS' ~ 1))

# Split the cART into their different groups
# Atripla
clinical_data_spec_var$ART1[clinical_data_spec_var$ART1 == "Atripla"] <- "Atripla_nrti"
clinical_data_spec_var$ART2[clinical_data_spec_var$ART1 == "Atripla_nrti"] <- "Atripla_nnrti"

# Genvoya
clinical_data_spec_var$ART1[clinical_data_spec_var$ART1 == "Genvoya"] <- "Genvoya_nrti"
clinical_data_spec_var$ART2[clinical_data_spec_var$ART1 == "Genvoya_nrti"] <- "Genvoya_pi"
clinical_data_spec_var$ART3[clinical_data_spec_var$ART1 == "Genvoya_pi"] <- "Genvoya_insti"

# Stribild
clinical_data_spec_var$ART1[clinical_data_spec_var$ART1 == "Stribild"] <- "Stribild_nrti"
clinical_data_spec_var$ART2[clinical_data_spec_var$ART1 == "Stribild_nrti"] <- "Stribild_pi"
clinical_data_spec_var$ART3[clinical_data_spec_var$ART1 == "Stribild_pi"] <- "Stribild_insti"

# Triumeq
clinical_data_spec_var$ART1[clinical_data_spec_var$ART1 == "Triumeq"] <- "Triumeq_nrti"
clinical_data_spec_var$ART2[clinical_data_spec_var$ART1 == "Triumeq_nrti"] <- "Triumeq_insti"


# List of unique names for all ART (curr + prev)
ART_all <- c(as_tibble(unique(clinical_data_spec_var$ART1))$value, as_tibble(unique(clinical_data_spec_var$ART2))$value, 
             as_tibble(unique(clinical_data_spec_var$ART3))$value, as_tibble(unique(clinical_data_spec_var$ART4))$value, 
             as_tibble(unique(clinical_data_spec_var$ART5))$value,
             as_tibble(unique(clinical_data_spec_var$ART1_prev))$value, as_tibble(unique(clinical_data_spec_var$ART2_prev))$value, 
             as_tibble(unique(clinical_data_spec_var$ART3_prev))$value, as_tibble(unique(clinical_data_spec_var$ART4_prev))$value, 
             as_tibble(unique(clinical_data_spec_var$ART5_prev))$value, as_tibble(unique(clinical_data_spec_var$ART6_prev))$value, 
             as_tibble(unique(clinical_data_spec_var$ART7_prev))$value, as_tibble(unique(clinical_data_spec_var$ART8_prev))$value, 
             as_tibble(unique(clinical_data_spec_var$ART9_prev))$value, as_tibble(unique(clinical_data_spec_var$ART10_prev))$value)
uniq_ART_all <- as_tibble(unique(ART_all))

# Based on the unique list of ART (43 terms), each term has been grouped into the following 6 groups
#The combined ART tablets: 'Atripla', 'Genvoya', 'Stribild', 'Triumeq' are split into the first 5 groups
NRTI            <- c('Combivir', 'Complera', 'Complera Complera', 'Descovy', 'Emtriva(emtricitabin,FTC)', 'Epivir(lamivudin,3TC)', 'Epzicom', 
                     'Kivexa', 'Lamivudin', 'Trizivir', 'Truvada', 'Viread(tenofovir,TDF)', 'Ziagen(abavacir,ABC)', 'Atripla_nrti', 'Genvoya_nrti', 'Stribild_nrti', 'Triumeq_nrti')
NNRTI           <- c('Edurant (rilpivirine, RPV)', 'Intelence (etravirine, ETR)', 'Invirase (saquinavir, SQV)', 'Stocrin', 
                     'Sustiva (Stocrin, efavirenz, EFV)', 'Viramune XR (nevirapine, NVP)', 'Atripla_nnrti')
PI              <- c('Amprenavir', 'Kaletra (Aluvia, lopinavir/ritonavir, LPV/r)', 'Norvir (ritonavir, RTV)', 'Norvir(ritonavir,RTV)','Prezcobix (Rezolsta, darunavir + cobicistat)', 
                     'Prezista (darunavir, DRV)', 'Reyataz (atazanavir, ATV)', 'Tybost(cobicstat)', 'Viracept (nelfinavir, NFV)', 'Genvoya_pi', 'Stribild_pi')
INSTI           <- c('Isentress (raltegravir)', 'Tivicay (dolutegravir)', 'Vitekta (elvitegravir)', 'Genvoya_insti', 'Stribild_insti', 'Triumeq_insti')
other_unknown   <- c('83', 'Fuzeon (enfuvirtide, ENF)', 'Selzentry (Celsentri, maraviroc)', 
                     'Retrovir(zidovudin,AZT)') # Only one lipo_art in "current ART cols", maybe a mistake - therefore in other group. Maraviroc is an entry inhibitor, however there are few of the, thus classified "other"

lipo_ART        <- c('Crixivan (indinavir, IDV)', 'Retrovir(zidovudin,AZT)', 'Videx(didanosine,ddl)', 'Videx(didanosine,ddl) Videx(didanosine)', 'Zerit(stavudin,d4T)')


# Creating new columns with encoded values for the ART groups and boolean column for thymidine exposure or not. 
druglist <- c('ART1', 'ART2', 'ART3', 'ART4', 'ART5')
clinical_data_spec_var_encoded <- clinical_data_spec_var %>% 
  mutate(NRTI = apply(.[, druglist], 1, function(r) any(r %in% NRTI)),
         NNRTI = apply(.[, druglist], 1, function(r) any(r %in% NNRTI)),
         PI = apply(.[, druglist], 1, function(r) any(r %in% PI)),
         INSTI = apply(.[, druglist], 1, function(r) any(r %in% INSTI)),
         Other_Unknown = apply(.[, druglist], 1, function(r) any(r %in% other_unknown)),
         thym_exposure = apply(.[, c('ART1_prev', 'ART2_prev', 'ART3_prev', 'ART4_prev', 'ART5_prev', 'ART6_prev', 'ART7_prev', 'ART8_prev', 'ART9_prev', 'ART10_prev')], 
                               1, function(r) any(r %in% lipo_ART)), # If a value from lipo_ART is in any of the 10 cols, then TRUE otherwise FALSE
         immunodeficiency = case_when(CD4_nadir > 200 | CDCAIDS == '0' ~ 0,
                                      CD4_nadir < 200 | CDCAIDS == '1' ~ 1)) %>% 
  select(ID_Condition, CONDITION, Condition_enc, AGE, GENDER, Gender_enc, Ethnic, immunodeficiency, thym_exposure, VAT, SAT,
         NRTI, NNRTI, PI, INSTI, Other_Unknown)


cols <- c('thym_exposure', 'NRTI', 'NNRTI', 'PI', 'INSTI', 'Other_Unknown')
clinical_data_spec_var_encoded[,cols] <- lapply(clinical_data_spec_var_encoded[,cols], as.numeric)






# Distribution of clinical data -------------------------------------------
#Compare GENDER -----------------------------------------------------------------
# Distribution of VAT stratified on condition
clinical_data_spec_var_encoded %>% 
  ggplot(mapping = aes(x = VAT, fill = CONDITION)) + 
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
  ggplot(mapping = aes(x = as.character(thym_exposure), fill = CONDITION)) + 
  geom_histogram(stat = "count", position = position_dodge()) +
  labs(title = 'Distribution of "Thymedine exposure" stratified on "Condition" divided in "Gender"', x = 'Thymidine exposure False(0) and True(1)') + 
  scale_fill_discrete(name='Condition',
                      breaks=c('HIV_NoMetS', 'HIV_MetS'),
                      labels=c('HIV without MetS', 'HIV with MetS')) +
  facet_grid(~ GENDER)

# Compare CONDITION -----------------------------------------------------------------
# Distribution of VAT stratified on thymidine exposure
clinical_data_spec_var_encoded %>% 
  ggplot(mapping = aes(x = VAT, fill = as.character(thym_exposure))) + 
  geom_density(alpha = 0.5) +
  labs(title = 'Distribution of "Waist-Hip ratio" stratified on thymidine exposure', x = 'Measure of "Waist to hip ratio"') +
  scale_fill_discrete(name='Thymidine Exposure',
                      breaks=c('0', '1'),
                      labels=c('False', 'True')) +
  facet_grid(~ CONDITION)

# Distribution of VAT stratified on age ranges
clinical_data_spec_var_encoded %>% 
  ggplot(mapping = aes(x = VAT, fill = as.character(age_groupings))) + 
  geom_density(alpha = 0.5) +
  labs(title = 'Distribution of "Waist-Hip ratio" stratified on age ranges', x = 'Measure of "Waist to hip ratio"') +
  scale_fill_discrete(name='Age groups',
                      breaks=c('1', '2', '3', '4'),
                      labels=c('40-49', '50-59', '60-69', '70-79')) +
  facet_grid(~ CONDITION)

# Distribution of VAT stratified on ART1 use
clinical_data_spec_var_encoded %>% 
  ggplot(mapping = aes(x = VAT, fill = as.character(NRTI))) + 
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

# Boxplot of VAT
clinical_data_spec_var_encoded %>% 
  ggplot(aes(GENDER, VAT, fill = CONDITION)) +
  geom_boxplot() 


# Boxplot of AGE
clinical_data_spec_var_encoded %>% 
  ggplot(aes(GENDER, AGE, fill = CONDITION)) +
  geom_boxplot() 


# Boxplot of "TAG54:4-FA20:3"
lipidomics_data %>% 
  filter(!Condition == "Ctrl") %>% 
  ggplot(aes(GENDER, `TAG54:4-FA20:3`, fill = Condition)) +
  geom_boxplot() 



# Drug plots
clinical_data_spec_var_encoded %>% 
  group_by(NRTI) %>% 
  count(.)  %>%  
  ggplot(aes(x = .$NRTI, .$n)) +
  geom_col()

clinical_data_spec_var_encoded %>% 
  group_by(NNRTI) %>% 
  count(.) %>% 
  ggplot(aes(x = .$NNRTI, .$n)) +
  geom_col()

clinical_data_spec_var_encoded %>% 
  group_by(PI) %>% 
  count(.)  %>%  
  ggplot(aes(x = .$PI, .$n)) +
  geom_col()

clinical_data_spec_var_encoded %>% 
  group_by(INSTI) %>% 
  count(.)  %>%  
  ggplot(aes(x = .$INSTI, .$n)) +
  geom_col()

clinical_data_spec_var_encoded %>% 
  group_by(Other_Unknown) %>% 
  count(.)  %>%  
  ggplot(aes(x = .$Other_Unknown, .$n)) +
  geom_col()









# Distribution of VAT stratified on ART1 use
clinical_data_spec_var_encoded %>% 
  ggplot(mapping = aes(x = VAT, fill = as.character(NRTI))) + 
  geom_histogram(alpha = 0.5) +
  labs(title = 'Distribution of "Waist-Hip ratio" stratified on ART1 use', x = 'Measure of "Waist to hip ratio"') 



# Summary statistics of clinical variables, difference between HIV_NoMetS vs. HIV_MetS -----------------------------------------------------------------
# Age
# Summary statistics
by(clinical_data_spec_var_encoded$AGE, clinical_data_spec_var_encoded$CONDITION, summary)
# Density plot, identifying the distribution
clinical_data_spec_var_encoded %>% 
  ggplot(mapping = aes(x = AGE)) + 
  geom_density(alpha = 0.5)
# Left skewed distribution --> Mann Whitney
wilcox.test(AGE ~ CONDITION, data = clinical_data_spec_var_encoded)

# Gender
# Summary statistics
by(clinical_data_spec_var_encoded$GENDER, clinical_data_spec_var_encoded$CONDITION, summary)
# Chi-squared test
chisq.test(clinical_data_spec_var_encoded$GENDER, clinical_data_spec_var_encoded$CONDITION)

# Ethnic
# Summary statistics
by(as.character(clinical_data_spec_var_encoded$Ethnic), clinical_data_spec_var_encoded$CONDITION, summary)
# Chi-squared test
chisq.test(as.character(clinical_data_spec_var_encoded$Ethnic), clinical_data_spec_var_encoded$CONDITION)

# CD4
# Density plot, identifying the distribution
clinical_data_spec_var_encoded %>% 
  ggplot(mapping = aes(x = CD4)) + 
  geom_density(alpha = 0.5) 
#Normality test: The data is normal if the p-value > 0.05. So the CD4 variable is NOT normally distributed.
ks.test(x = clinical_data_spec_var_encoded$CD4, y = "pnorm", mean = mean(clinical_data_spec_var_encoded$CD4), sd = sd(clinical_data_spec_var_encoded$CD4))$p.value #p-value = 1.014642e-08
# Left skewed distribution --> Mann Whitney
wilcox.test(CD4 ~ CONDITION, data = clinical_data_spec_var_encoded)
# Summary statistics
by(clinical_data_spec_var_encoded$CD4, clinical_data_spec_var_encoded$CONDITION, summary)

# VAT
clinical_data_spec_var_encoded %>% 
  ggplot(mapping = aes(x = VAT)) + 
  geom_density(alpha = 0.5)
# Normal distribution --> T-test
t.test(VAT ~ CONDITION, data = clinical_data_spec_var_encoded)
# Summary statistics
by(clinical_data_spec_var_encoded$VAT, clinical_data_spec_var_encoded$CONDITION, summary)


# SAT
clinical_data_spec_var_encoded %>% 
  ggplot(mapping = aes(x = SAT)) + 
  geom_density(alpha = 0.5)
# Normal distribution --> T-test
t.test(SAT ~ CONDITION, data = clinical_data_spec_var_encoded)
# Summary statistics
by(clinical_data_spec_var_encoded$SAT, clinical_data_spec_var_encoded$CONDITION, summary)


# Immunodeficency
# Summary statistics
by(as.character(clinical_data_spec_var_encoded$immunodeficiency), clinical_data_spec_var_encoded$CONDITION, summary)
# Chi-squared test
chisq.test(as.character(clinical_data_spec_var_encoded$immunodeficiency), clinical_data_spec_var_encoded$CONDITION)

# Thymidine exposure
# Summary statistics
by(as.character(clinical_data_spec_var_encoded$thym_exposure), clinical_data_spec_var_encoded$CONDITION, summary)
# Chi-squared test
chisq.test(as.character(clinical_data_spec_var_encoded$thym_exposure), clinical_data_spec_var_encoded$CONDITION)

# NRTI
# Summary statistics
by(as.character(clinical_data_spec_var_encoded$NRTI), clinical_data_spec_var_encoded$CONDITION, summary)
# Fisher t-test --> small sample size for INSTI group
fisher.test(as.character(clinical_data_spec_var_encoded$NRTI), clinical_data_spec_var_encoded$CONDITION)

# NNRTI
# Summary statistics
by(as.character(clinical_data_spec_var_encoded$NNRTI), clinical_data_spec_var_encoded$CONDITION, summary)
# Fisher t-test --> small sample size
fisher.test(as.character(clinical_data_spec_var_encoded$NNRTI), clinical_data_spec_var_encoded$CONDITION)

# PI
# Summary statistics
by(as.character(clinical_data_spec_var_encoded$PI), clinical_data_spec_var_encoded$CONDITION, summary)
# Fisher t-test --> small sample size
fisher.test(as.character(clinical_data_spec_var_encoded$PI), clinical_data_spec_var_encoded$CONDITION)

# INSTI
# Summary statistics
by(as.character(clinical_data_spec_var_encoded$INSTI), clinical_data_spec_var_encoded$CONDITION, summary)
# Fisher t-test --> small sample size
fisher.test(as.character(clinical_data_spec_var_encoded$INSTI), clinical_data_spec_var_encoded$CONDITION)


# Other_Unknown 
# Summary statistics
by(as.character(clinical_data_spec_var_encoded$Other_Unknown), clinical_data_spec_var_encoded$CONDITION, summary)
# Fisher t-test --> small sample size
fisher.test(as.character(clinical_data_spec_var_encoded$Other_Unknown), clinical_data_spec_var_encoded$CONDITION)







# Utables -----------------------------------------------------------------
# Change class from numeric to chracther --> should be moved up in the script....
cols = c('Ethnic', 'immunodeficiency', 'thym_exposure','NRTI', 'NNRTI', 'PI', 'INSTI', 'Other_Unknown');    
clinical_data_spec_var_encoded[,cols] = apply(clinical_data_spec_var_encoded[,cols], 2, function(x) as.character(as.numeric(x)));

# Check class type
sapply(clinical_data_spec_var_encoded, class)

# Use utable to characterise data
univariateTable(CONDITION ~ 
                  GENDER+AGE+Ethnic+immunodeficiency+thym_exposure+VAT+SAT+
                  NRTI+NNRTI+PI+INSTI+Other_Unknown, 
                data = clinical_data_spec_var_encoded)



# Regression models -------------------------------------------------------
# Merge all lipids and significant metabolites conc.
lipidomics_met_data <- merge(lipidomics_data, metabolite_sign, by = "ID_Condition") %>% t()
colnames(lipidomics_met_data) <- lipidomics_data$ID_Condition

# Include only lipids/metabolites from positive weigthed network
common_cols <- intersect(network_table$Biochemicals, rownames(lipidomics_met_data))

# Find outersect
#a <- colnames(network_table_t)
#b <- colnames(lipidomics_met_data[, 4:ncol(lipidomics_met_data)])
#setdiff(union(a,b), intersect(a,b))

# 
network_table_2 <- as.data.frame(network_table) %>% filter(Biochemicals %in% common_cols)
lipidomics_met_data_2 <- as.data.frame(lipidomics_met_data) %>% 
  mutate(Biochemicals = rownames(lipidomics_met_data)) %>% 
  filter(Biochemicals %in% common_cols) %>% 
  select(Biochemicals, everything())

combined_table <- merge(network_table_2, lipidomics_met_data_2, by = "Biochemicals") %>% 
  select(everything(), -contains("Ctrl"))

# From factor to numeric
samples <- colnames(combined_table[9:ncol(combined_table)])
combined_table <- as_tibble(combined_table)
combined_table[,samples] = apply(combined_table[,samples], 2, function(x) as.numeric(as.character(x)));

# Function calculating z-score on each biochemical cell: 
# The z-score is a measure that shows how much away (below or above) of the mean is a specific value (individual) in a given dataset.
zscore <- function(x) {
  (x-mean(x))/sd(x)
}

# Compute the z-score for each lipid
piv_ <- combined_table %>% 
  select(everything(), -c(degree, logFC_limma, adj_pval_limma, Group, Super_pathway, community, sign_lip_met)) %>% 
  pivot_longer(-Biochemicals) %>% 
  pivot_wider(names_from=Biochemicals, values_from=value)

combined_zscore <- cbind(piv_[1],lapply(piv_[2:ncol(piv_)], zscore))

combined_zscore_2 <- combined_zscore %>% 
  rename(Biochemicals = name) %>% 
  pivot_longer(-Biochemicals) %>% 
  pivot_wider(names_from=Biochemicals, values_from=value) 

combined_zscore <- cbind(combined_table[1:8], combined_zscore_2) %>% select(everything(), -name)

# Compute the average of the z-score for each community
community_score <- combined_zscore %>% 
  group_by(community) %>% 
  summarise_at(samples, mean, na.rm = TRUE) %>% 
  pivot_longer(-community) %>% 
  pivot_wider(names_from=community, values_from=value) 


# Encode clinical variables
cols = c('immunodeficiency', 'thym_exposure');    
clinical_data_spec_var_encoded[,cols] = apply(clinical_data_spec_var_encoded[,cols], 2, function(x) as.numeric(as.character(x)));

clinical_data_spec_var_encoded_2 <- clinical_data_spec_var_encoded[,c("ID_Condition", "Condition_enc","AGE", "Ethnic","Gender_enc",
                                                                      "immunodeficiency","thym_exposure", "VAT", "SAT", 
                                                                      "NRTI", "NNRTI", "PI", "INSTI", "Other_Unknown")]

# Make a table including clincal variables and community variables
clinical_community_data <- cbind(clinical_data_spec_var_encoded_2, community_score) %>% 
  select(everything(), -name)

# Calculate the average of the z-score for each community
combined_zscore_2 <- combined_zscore %>%
  select(everything(), -c(degree, logFC_limma, adj_pval_limma, Group, Super_pathway, community, sign_lip_met)) %>% 
  pivot_longer(-Biochemicals) %>% 
  pivot_wider(names_from=Biochemicals, values_from=value) 

clinical_community_values <- cbind(clinical_data_spec_var_encoded_2, community_score) %>% select(everything(), -name)

regression_data <- cbind(clinical_community_values, combined_zscore_2) %>% select(everything(), -name)



# For VAT
lmmodel1 = lm(VAT~AGE+Gender_enc+Condition_enc+NRTI+NNRTI+`TAG54:4-FA20:3`, data = regression_data) #Create the linear regression
summary(lmmodel1)
summary(lmmodel1)$coefficient
##publish(lmmodel1)
##summary(regressionTable(lmmodel1)) #Review the results

plot(lmmodel1, pch = 16, col = "blue") #Plot the results

qqnorm(lmmodel1$residuals)
qqline(lmmodel1$residuals)





# Model 1 -----------------------------------------------------------------
# Clin_x = age + sex + condition + lipid_x
# Check how each lipid respond to each variable, this is supported by heatmap of principal component correlated with clinical variables

# Function creating a coefficient table of lipids and metabolites
coeff_table <- function(multi_lms) {
  # summary
  multi_lms_summaries <- lapply(multi_lms, summary)
  multi_lms_coef <- lapply(multi_lms_summaries, function(x) x$coefficients)
  
  # Loop extracting only the lipid values from the regression model
  biochemical_l <- list()
  for(i in multi_lms_coef) {
    new_element <- list(i[5,])                       
    biochemical_l[[length(biochemical_l) + 1]] <- new_element
  }
  
  model_biochemicals <- as.data.frame(biochemical_l) %>% 
    t(.) %>% cbind(combined_zscore$Biochemicals, .) %>%  
    as_tibble(.) %>% 
    rename(Biochemicals = V1)
  
  model <- merge(model_biochemicals, network_table, by = "Biochemicals") %>% 
    as_tibble(.) %>% 
    select(Biochemicals, Estimate, `Std. Error`, `t value`, `Pr(>|t|)`, community, sign_lip_met) 
  return(model)
}

biochemical_start <- 18

# Run multiple regressions of all lipids/metabolites, adjusting for age, sex, condition
# VAT
multi_lms_VAT <- lapply(colnames(regression_data[,biochemical_start:ncol(regression_data)]), function(x) lm(VAT~AGE+Gender_enc+Condition_enc+regression_data[,x], data = regression_data))
model_VAT <- coeff_table(multi_lms_VAT) %>% 
  mutate(pvalue_adjusted = round(p.adjust(.$`Pr(>|t|)`, "fdr"),5))

# SAT
multi_lms_SAT <- lapply(colnames(regression_data[,biochemical_start:ncol(regression_data)]), function(x) lm(SAT~AGE+Gender_enc+Condition_enc+regression_data[,x], data = regression_data))
model_SAT <- coeff_table(multi_lms_SAT) %>% 
  mutate(pvalue_adjusted = round(p.adjust(.$`Pr(>|t|)`, "fdr"),5))

# immunodecifiency
multi_lms_immunodecifiency <- lapply(colnames(regression_data[,biochemical_start:ncol(regression_data)]), function(x) lm(immunodeficiency~AGE+Gender_enc+Condition_enc+regression_data[,x], data = regression_data))
model_immunodecifiency <- coeff_table(multi_lms_immunodecifiency)  %>% 
  mutate(pvalue_adjusted = p.adjust(.$`Pr(>|t|)`, "fdr"))

# Ethnic
multi_lms_Ethnic <- lapply(colnames(regression_data[,biochemical_start:ncol(regression_data)]), function(x) lm(Ethnic~AGE+Gender_enc+Condition_enc+regression_data[,x], data = regression_data))
model_Ethnic <- coeff_table(multi_lms_Ethnic)  %>% 
  mutate(pvalue_adjusted = p.adjust(.$`Pr(>|t|)`, "fdr"))

# thymidine exposure
multi_lms_thym_exposure <- lapply(colnames(regression_data[,biochemical_start:ncol(regression_data)]), function(x) lm(thym_exposure~AGE+Gender_enc+Condition_enc+regression_data[,x], data = regression_data))
model_thym_exposure <- coeff_table(multi_lms_thym_exposure)  %>% 
  mutate(pvalue_adjusted = round(p.adjust(.$`Pr(>|t|)`, "fdr"),5))

# NRTI 
multi_lms_NRTI <- lapply(colnames(regression_data[,biochemical_start:ncol(regression_data)]), function(x) lm(NRTI~AGE+Gender_enc+Condition_enc+regression_data[,x], data = regression_data))
model_NRTI <- coeff_table(multi_lms_NRTI)  %>% 
  mutate(pvalue_adjusted = round(p.adjust(.$`Pr(>|t|)`, "fdr"),5))

# NNRTI 
multi_lms_NNRTI <- lapply(colnames(regression_data[,biochemical_start:ncol(regression_data)]), function(x) lm(NNRTI~AGE+Gender_enc+Condition_enc+regression_data[,x], data = regression_data))
model_NNRTI <- coeff_table(multi_lms_NNRTI)  %>% 
  mutate(pvalue_adjusted = round(p.adjust(.$`Pr(>|t|)`, "fdr"),5))

# PI 
multi_lms_PI <- lapply(colnames(regression_data[,biochemical_start:ncol(regression_data)]), function(x) lm(PI~AGE+Gender_enc+Condition_enc+regression_data[,x], data = regression_data))
model_PI <- coeff_table(multi_lms_PI)  %>% 
  mutate(pvalue_adjusted = round(p.adjust(.$`Pr(>|t|)`, "fdr"),5))

# INSTI 
multi_lms_INSTI <- lapply(colnames(regression_data[,biochemical_start:ncol(regression_data)]), function(x) lm(INSTI~AGE+Gender_enc+Condition_enc+regression_data[,x], data = regression_data))
model_INSTI <- coeff_table(multi_lms_INSTI)  %>% 
  mutate(pvalue_adjusted = round(p.adjust(.$`Pr(>|t|)`, "fdr"),5))


# Model 2 -----------------------------------------------------------------
# Clin_x = z-score
# Not adjusting for any variable

# Function creating a coefficient table of comminities
comm_coeff_table <- function(multi_lms) {
  # summary
  multi_lms_summaries <- lapply(multi_lms, summary)
  multi_lms_coef <- lapply(multi_lms_summaries, function(x) x$coefficients)
  
  # Loop extracting only the community values from the regression model
  comm_l <- list()
  for(i in multi_lms_coef) {
    new_element <- list(i[2,])                       
    comm_l[[length(comm_l) + 1]] <- new_element
  }
  
  model_community <- as.data.frame(comm_l) %>% 
    t(.) %>% cbind(unique(combined_zscore$community), .) %>%  
    as_tibble(.)  %>% 
    rename(community = V1) #%>% 
  #filter(community %in% c('c1', 'c2', 'c3')) 
  return(model_community)
}

comm_cols = c('c1', 'c2', 'c3') #, 'c4', 'c5', 'c6', 'c7')

# Age
community_lms_Age <- lapply(colnames(regression_data[,comm_cols]), function(x) lm(AGE~regression_data[,x], data = regression_data))
model_comm_Age <- comm_coeff_table(community_lms_Age) %>% 
  mutate(pvalue_adjusted = round(p.adjust(.$`Pr(>|t|)`, "fdr"),5))

# Condition
community_lms_condition <- lapply(colnames(regression_data[,comm_cols]), function(x) lm(Condition_enc~regression_data[,x], data = regression_data))
model_comm_condition <- comm_coeff_table(community_lms_condition) %>% 
  mutate(pvalue_adjusted = round(p.adjust(.$`Pr(>|t|)`, "fdr"),5))

# Sex
community_lms_sex <- lapply(colnames(regression_data[,comm_cols]), function(x) lm(Gender_enc~regression_data[,x], data = regression_data))
model_comm_sex <- comm_coeff_table(community_lms_sex) %>% 
  mutate(pvalue_adjusted = round(p.adjust(.$`Pr(>|t|)`, "fdr"),5))

# Ethnic
community_lms_ethnic <- lapply(colnames(regression_data[,comm_cols]), function(x) lm(Ethnic~regression_data[,x], data = regression_data))
model_comm_ethnic <- comm_coeff_table(community_lms_ethnic) %>% 
  mutate(pvalue_adjusted = round(p.adjust(.$`Pr(>|t|)`, "fdr"),5))

# VAT
community_lms_VAT <- lapply(colnames(regression_data[,comm_cols]), function(x) lm(VAT~regression_data[,x], data = regression_data))
model_comm_VAT <- comm_coeff_table(community_lms_VAT) %>% 
  mutate(pvalue_adjusted = round(p.adjust(.$`Pr(>|t|)`, "fdr"),5))

# SAT
community_lms_SAT <- lapply(colnames(regression_data[,comm_cols]), function(x) lm(SAT~regression_data[,x], data = regression_data))
model_comm_SAT <- comm_coeff_table(community_lms_SAT) %>% 
  mutate(pvalue_adjusted = round(p.adjust(.$`Pr(>|t|)`, "fdr"),5))

# immunodeficiency
community_lms_immunodeficiency <- lapply(colnames(regression_data[,comm_cols]), function(x) lm(immunodeficiency~regression_data[,x], data = regression_data))
model_comm_immunodeficiency <- comm_coeff_table(community_lms_immunodeficiency) %>% 
  mutate(pvalue_adjusted = round(p.adjust(.$`Pr(>|t|)`, "fdr"),5))

# thymidine exposure
community_lms_thym_exposure <- lapply(colnames(regression_data[,comm_cols]), function(x) lm(thym_exposure~regression_data[,x], data = regression_data))
model_comm_thym_exposure <- comm_coeff_table(community_lms_thym_exposure) %>% 
  mutate(pvalue_adjusted = round(p.adjust(.$`Pr(>|t|)`, "fdr"),5))


# NRTI 
community_lms_NRTI <- lapply(colnames(regression_data[,comm_cols]), function(x) lm(NRTI~regression_data[,x], data = regression_data))
model_comm_NRTI <- comm_coeff_table(community_lms_NRTI) %>% 
  mutate(pvalue_adjusted = round(p.adjust(.$`Pr(>|t|)`, "fdr"),5))

# NNRTI 
community_lms_NNRTI <- lapply(colnames(regression_data[,comm_cols]), function(x) lm(NNRTI~regression_data[,x], data = regression_data))
model_comm_NNRTI <- comm_coeff_table(community_lms_NNRTI)   %>% 
  mutate(pvalue_adjusted = round(p.adjust(.$`Pr(>|t|)`, "fdr"),5))

# PI 
community_lms_PI <- lapply(colnames(regression_data[,comm_cols]), function(x) lm(PI~regression_data[,x], data = regression_data))
model_comm_PI <- comm_coeff_table(community_lms_PI)   %>% 
  mutate(pvalue_adjusted = round(p.adjust(.$`Pr(>|t|)`, "fdr"),5))

# INSTI 
community_lms_INSTI <- lapply(colnames(regression_data[,comm_cols]), function(x) lm(INSTI~regression_data[,x], data = regression_data))
model_comm_INSTI <- comm_coeff_table(community_lms_INSTI)   %>% 
  mutate(pvalue_adjusted = round(p.adjust(.$`Pr(>|t|)`, "fdr"),5))



# Heatmap of PC and clinical variables ------------------------------------

library(Biobase)
library(PCAtools)

rownames(clinical_data_spec_var_encoded) <- clinical_data_spec_var_encoded$ID_Condition
lipidomics_met_data_3 <- cbind(lipidomics_met_data_2[1], apply(lipidomics_met_data_2[,2:ncol(lipidomics_met_data_2)], 2, function(x) as.numeric(as.factor(x)))) %>% select(rownames(clinical_data_spec_var_encoded))

# PCA object 
p <- pca(as.data.frame(lipidomics_met_data_3), metadata = as.data.frame(clinical_data_spec_var_encoded), removeVar = 0.1)

# Scree plot, including 90% of the explained variance
screeplot(p, components = getComponents(p)[0:30], axisLabSize = 18, titleLabSize = 22)

eigencorplot(p, metavars = c("Condition_enc", "Ethnic", "AGE","Gender_enc", "VAT", "SAT", 
                             "immunodeficiency","thym_exposure",
                             "NRTI", "NNRTI", "PI", "INSTI", "Other_Unknown"))



biplot(p, colby = 'SAT',
       hline = 0, vline = 0,
       legendPosition = 'right')

#OBS naar du kommer tilbage skal du lave PCA og begynde paa p, hvor du siger at dette community associerede med disse kliniske variable osv.



# Optimal number of PCs, illustrated by scree plot
horn <- parallelPCA(as.data.frame(lipidomics_met_data_3))
elbow <- findElbowPoint(p$variance)

screeplot(p,components = getComponents(p, 1:20),
          vline = c(horn$n, elbow)) +
  geom_label(aes(x = horn$n + 1, y = 50,
                 label = 'Horn\'s', vjust = -1, size = 8)) +
  geom_label(aes(x = elbow + 1, y = 50,
                 label = 'Elbow method', vjust = -1, size = 8))

# Check how many PC covers X % of the variance
which(cumsum(p$variance) > 80)[1] # n = 13


# Correlate the principal components back to the clinical data
eigencorplot(p,components = getComponents(p)[1:13],
             metavars = c("Condition_enc","AGE", "Ethnic", "Gender_enc","VAT","SAT",
                          "immunodeficiency","thym_exposure",
                          "NRTI", "NNRTI", "PI", "INSTI", "Other_Unknown"),
             cexCorval = 0.7,
             colCorval = 'white',
             fontCorval = 2,
             posLab = 'bottomleft',
             rotLabX = 45,
             posColKey = 'top',
             cexLabColKey = 1.5,
             scale = TRUE,
             main = 'Correlation between PC[1-13] (80% explained var) and the clinical variables',
             colFrame = 'white',
             plotRsquared = FALSE)
