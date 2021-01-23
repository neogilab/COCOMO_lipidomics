##########################################################
# Script for regression analysis
##########################################################

# Clear workspace
rm(list = ls())

# Load packages
library(tidyverse)
library(patchwork)
library(PCAtools)

# Load data 
lipidomics_data <- read_csv("data/02_lipidomics_data_tidy.csv")
metabolite_sign <- read_csv("data/02_metabolomics_sign_metabolites.csv")
clinical_data <- read_csv("data/02_clinical_data_tidy.csv")
network_table <- read_csv("data/08_cytoscape_pos_table_FDR.csv") 



# Data preparation -------------------------------------------------------
# Merge all lipids and significant metabolites conc.
lipidomics_met_data <- merge(lipidomics_data, metabolite_sign, by = "ID_Condition") %>% t()
colnames(lipidomics_met_data) <- lipidomics_data$ID_Condition

# Include only lipids/metabolites from positive weigthed network
common_cols <- intersect(network_table$Biochemicals, rownames(lipidomics_met_data))

# Make table with lipids and metabolites from the pos network
###network_table_2 <- as.data.frame(network_table) %>% filter(Biochemicals %in% common_cols)
lipidomics_met_data_2 <- as.data.frame(lipidomics_met_data) %>% 
  mutate(Biochemicals = rownames(lipidomics_met_data)) %>% 
  filter(Biochemicals %in% common_cols) %>% 
  select(Biochemicals, everything())

combined_table <- merge(network_table, lipidomics_met_data_2, by = "Biochemicals") %>% 
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
clinical_data[,cols] = apply(clinical_data[,cols], 2, function(x) as.numeric(as.character(x)));

clinical_data_2 <- clinical_data[,c("ID_Condition", "Condition_enc","AGE", "Ethnic","Gender_enc",
                                    "immunodeficiency", "thym_exposure", "VAT", "SAT", 
                                    "NRTI", "NNRTI", "PI", "INSTI", "Other_Unknown")]

# Make a table including clincal variables and community variables
clinical_community_data <- cbind(clinical_data_2, community_score) %>% 
  select(everything(), -name)

# Calculate the average of the z-score for each community
combined_zscore_2 <- combined_zscore %>%
  select(everything(), -c(degree, logFC_limma, adj_pval_limma, Group, Super_pathway, community, sign_lip_met)) %>% 
  pivot_longer(-Biochemicals) %>% 
  pivot_wider(names_from=Biochemicals, values_from=value) 

clinical_community_values <- cbind(clinical_data_2, community_score) %>% select(everything(), -name)
regression_data <- cbind(clinical_community_values, combined_zscore_2) %>% select(everything(), -name)



# Heatmap of PC vs. clinical variables ------------------------------------
# Check the trend of how the whole lipid/metabolite network correlates with the clinical variables
clinical_data_2 <- clinical_data %>% 
  select(MetS = Condition_enc, 
         Age = AGE, 
         Sex = Gender_enc,
         Ethnicity = Ethnic,
         Immunodeficiency = immunodeficiency,
         `Exposure to preART` = thym_exposure, 
         VAT, 
         SAT, 
         `ART_NRTI` = NRTI, 
         `ART_NNRTI` = NNRTI, 
         `ART_PI` = PI, 
         `ART_INSTI` = INSTI,
         `ART_other/unknown` = Other_Unknown
  )
rownames(clinical_data_2) <- clinical_data$ID_Condition

# Changing from factor to numeric
lipidomics_met_data_3 <- lipidomics_met_data_2 %>% select(rownames(clinical_data_2)) 
#indx <- sapply(lipidomics_met_data_3, is.factor)
indx <- colnames(lipidomics_met_data_3)
lipidomics_met_data_3[,indx] = apply(lipidomics_met_data_3[,indx], 2, function(x) as.numeric(as.character(x)));

# PCA object 
p <- pca(as.data.frame(lipidomics_met_data_3), metadata = as.data.frame(clinical_data_2), removeVar = 0.1)

# Scree plot, including 90% of the explained variance - Optimal number of PCs, illustrated by scree plot
horn <- parallelPCA(as.data.frame(lipidomics_met_data_3))
elbow <- findElbowPoint(p$variance)

png("results/10_1_clinical_vs_omics_PC_screeplot.png",#,
    width = 6*300,        # 5 x 300 pixels
    height = 6*300,
    res = 200)            # font size)    
screeplot(p,components = getComponents(p, 1:10))
dev.off()

# Check how many PC covers X % of the variance
which(cumsum(p$variance) > 90)[1] # n = 6

# Correlation between PC[1-6] (90% explained var) and the clinical variables
png("results/10_2_clinical_vs_omics_PC_heatmap.png",#,
    width = 6*300,        # 5 x 300 pixels
    height = 6*300,
    res = 200,            # 300 pixels per inch
    pointsize = 150)      # font size)    
eigencorplot(p,components = getComponents(p)[1:6],
             metavars = c("MetS","Age","Sex","Ethnicity","Immunodeficiency","Exposure to preART","VAT" , 
                          "SAT","ART_NRTI","ART_NNRTI","ART_PI","ART_INSTI","ART_other/unknown"),
             cexCorval = 0.7,
             colCorval = 'white',
             fontCorval = 2,
             posLab = 'bottomleft',
             rotLabX = 45,
             posColKey = 'top',
             cexLabColKey = 1.5,
             scale = TRUE,
             colFrame = 'white',
             plotRsquared = FALSE,
             corMultipleTestCorrection = "BH")
dev.off()



# Regression model 1 -----------------------------------------------------------------
# Check the trend of how each community correlates with the clinical variables, not adjusting for any variable
# Clin_x = z-score

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

comm_cols = c('c1', 'c2', 'c3') 

# Age
community_lms_Age <- lapply(colnames(regression_data[,comm_cols]), function(x) lm(AGE~regression_data[,x], data = regression_data))
model_comm_Age <- comm_coeff_table(community_lms_Age) %>% 
  mutate(clinical_var = "age") %>% 
  select(community, clinical_var, everything())

# MetS
community_lms_condition <- lapply(colnames(regression_data[,comm_cols]), function(x) lm(Condition_enc~regression_data[,x], data = regression_data))
model_comm_condition <- comm_coeff_table(community_lms_condition) %>% 
  mutate(clinical_var = "condition") %>% 
  select(community, clinical_var, everything())

# Sex
community_lms_sex <- lapply(colnames(regression_data[,comm_cols]), function(x) lm(Gender_enc~regression_data[,x], data = regression_data))
model_comm_sex <- comm_coeff_table(community_lms_sex) %>% 
  mutate(clinical_var = "sex") %>% 
  select(community, clinical_var, everything())

# Ethnicity
community_lms_Ethnicity <- lapply(colnames(regression_data[,comm_cols]), function(x) lm(Ethnic~regression_data[,x], data = regression_data))
model_comm_Ethnicity <- comm_coeff_table(community_lms_Ethnicity) %>% 
  mutate(clinical_var = "ethnicity") %>% 
  select(community, clinical_var, everything())

# VAT
community_lms_VAT <- lapply(colnames(regression_data[,comm_cols]), function(x) lm(VAT~regression_data[,x], data = regression_data))
model_comm_VAT <- comm_coeff_table(community_lms_VAT) %>% 
  mutate(clinical_var = "VAT") %>% 
  select(community, clinical_var, everything())

# SAT
community_lms_SAT <- lapply(colnames(regression_data[,comm_cols]), function(x) lm(SAT~regression_data[,x], data = regression_data))
model_comm_SAT <- comm_coeff_table(community_lms_SAT) %>% 
  mutate(clinical_var = "SAT") %>% 
  select(community, clinical_var, everything())

# Immunodeficiency
community_lms_Immunodeficiency <- lapply(colnames(regression_data[,comm_cols]), function(x) lm(immunodeficiency~regression_data[,x], data = regression_data))
model_comm_Immunodeficiency <- comm_coeff_table(community_lms_Immunodeficiency) %>% 
  mutate(clinical_var = "immunodeficiency") %>% 
  select(community, clinical_var, everything())

# Exposure to preART
community_lms_thym_exposure <- lapply(colnames(regression_data[,comm_cols]), function(x) lm(thym_exposure~regression_data[,x], data = regression_data))
model_comm_thym_exposure <- comm_coeff_table(community_lms_thym_exposure) %>% 
  mutate(clinical_var = "exposure to preART") %>% 
  select(community, clinical_var, everything())

# NRTI 
community_lms_NRTI <- lapply(colnames(regression_data[,comm_cols]), function(x) lm(NRTI~regression_data[,x], data = regression_data))
model_comm_NRTI <- comm_coeff_table(community_lms_NRTI) %>% 
  mutate(clinical_var = "ART_NRTI") %>% 
  select(community, clinical_var, everything())

# NNRTI 
community_lms_NNRTI <- lapply(colnames(regression_data[,comm_cols]), function(x) lm(NNRTI~regression_data[,x], data = regression_data))
model_comm_NNRTI <- comm_coeff_table(community_lms_NNRTI)   %>% 
  mutate(clinical_var = "ART_NNRTI") %>% 
  select(community, clinical_var, everything())

# PI 
community_lms_PI <- lapply(colnames(regression_data[,comm_cols]), function(x) lm(PI~regression_data[,x], data = regression_data))
model_comm_PI <- comm_coeff_table(community_lms_PI)   %>% 
  mutate(clinical_var = "ART_PI") %>% 
  select(community, clinical_var, everything())

# INSTI 
community_lms_INSTI <- lapply(colnames(regression_data[,comm_cols]), function(x) lm(INSTI~regression_data[,x], data = regression_data))
model_comm_INSTI <- comm_coeff_table(community_lms_INSTI)   %>% 
  mutate(clinical_var = "ART_INSTI") %>% 
  select(community, clinical_var, everything())

# ART_other/unknown 
community_lms_Other_Unknown <- lapply(colnames(regression_data[,comm_cols]), function(x) lm(Other_Unknown~regression_data[,x], data = regression_data))
model_comm_Other_Unknown<- comm_coeff_table(community_lms_Other_Unknown)   %>% 
  mutate(clinical_var = "ART_other_unknown") %>% 
  select(community, clinical_var, everything())

# Adjust with FDR
all_stats_community_clinical <- rbind(model_comm_Age, model_comm_condition, model_comm_Ethnicity, model_comm_Immunodeficiency, 
                                      model_comm_INSTI, model_comm_NNRTI, model_comm_NRTI, model_comm_PI, model_comm_SAT, 
                                      model_comm_sex, model_comm_thym_exposure, model_comm_VAT, model_comm_Other_Unknown) %>% 
  mutate(pvalue_adjusted = round(p.adjust(.$`Pr(>|t|)`, "fdr"),5))

# From character to numeric
cols = c("Estimate","Std. Error", "t value", "Pr(>|t|)", "pvalue_adjusted")    
all_stats_community_clinical[,cols] = apply(all_stats_community_clinical[,cols], 2, function(x) as.numeric(as.character(x)));

model1_all_stats_community_clinical <- cbind(all_stats_community_clinical[, c("community", "clinical_var")], 
                                             round(all_stats_community_clinical[,c("Estimate","Std. Error", "t value")],2), 
                                             round(all_stats_community_clinical[, c("Pr(>|t|)", "pvalue_adjusted")], 6)) %>% 
  write_csv("data/model1_community_vs_clinical.csv")

# Regression model 2 -----------------------------------------------------------------
# Check the trend of how each lipid and key metabolite correlates with the clinical variables
# Clin_x = age + sex + condition + lipid_x

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
  mutate(clinical_var = "VAT") %>% 
  select(community, clinical_var, everything())

# SAT
multi_lms_SAT <- lapply(colnames(regression_data[,biochemical_start:ncol(regression_data)]), function(x) lm(SAT~AGE+Gender_enc+Condition_enc+regression_data[,x], data = regression_data))
model_SAT <- coeff_table(multi_lms_SAT)  %>%   
  mutate(clinical_var = "SAT") %>% 
  select(community, clinical_var, everything())

# immunodecifiency
multi_lms_immunodecifiency <- lapply(colnames(regression_data[,biochemical_start:ncol(regression_data)]), function(x) lm(immunodeficiency~AGE+Gender_enc+Condition_enc+regression_data[,x], data = regression_data))
model_immunodecifiency <- coeff_table(multi_lms_immunodecifiency)   %>%   
  mutate(clinical_var = "Immunodeficiency") %>% 
  select(community, clinical_var, everything())

# Ethnicity
multi_lms_Ethnicity <- lapply(colnames(regression_data[,biochemical_start:ncol(regression_data)]), function(x) lm(Ethnic~AGE+Gender_enc+Condition_enc+regression_data[,x], data = regression_data))
model_Ethnicity <- coeff_table(multi_lms_Ethnicity)  %>%   
  mutate(clinical_var = "Ethnicity") %>% 
  select(community, clinical_var, everything())

# thymidine exposure
multi_lms_thym_exposure <- lapply(colnames(regression_data[,biochemical_start:ncol(regression_data)]), function(x) lm(thym_exposure~AGE+Gender_enc+Condition_enc+regression_data[,x], data = regression_data))
model_thym_exposure <- coeff_table(multi_lms_thym_exposure)  %>%   
  mutate(clinical_var = "Exposure to preART") %>% 
  select(community, clinical_var, everything())

# NRTI 
multi_lms_NRTI <- lapply(colnames(regression_data[,biochemical_start:ncol(regression_data)]), function(x) lm(NRTI~AGE+Gender_enc+Condition_enc+regression_data[,x], data = regression_data))
model_NRTI <- coeff_table(multi_lms_NRTI)   %>%   
  mutate(clinical_var = "ART_NRTI") %>% 
  select(community, clinical_var, everything())

# NNRTI 
multi_lms_NNRTI <- lapply(colnames(regression_data[,biochemical_start:ncol(regression_data)]), function(x) lm(NNRTI~AGE+Gender_enc+Condition_enc+regression_data[,x], data = regression_data))
model_NNRTI <- coeff_table(multi_lms_NNRTI)  %>%   
  mutate(clinical_var = "ART_NNRTI") %>% 
  select(community, clinical_var, everything())

# PI 
multi_lms_PI <- lapply(colnames(regression_data[,biochemical_start:ncol(regression_data)]), function(x) lm(PI~AGE+Gender_enc+Condition_enc+regression_data[,x], data = regression_data))
model_PI <- coeff_table(multi_lms_PI)  %>%   
  mutate(clinical_var = "ART_PI") %>% 
  select(community, clinical_var, everything())

# INSTI 
multi_lms_INSTI <- lapply(colnames(regression_data[,biochemical_start:ncol(regression_data)]), function(x) lm(INSTI~AGE+Gender_enc+Condition_enc+regression_data[,x], data = regression_data))
model_INSTI <- coeff_table(multi_lms_INSTI)   %>%   
  mutate(clinical_var = "ART_INSTI") %>% 
  select(community, clinical_var, everything())

# ART_other/unknown 
multi_lms_other_unknown <- lapply(colnames(regression_data[,biochemical_start:ncol(regression_data)]), function(x) lm(Other_Unknown~AGE+Gender_enc+Condition_enc+regression_data[,x], data = regression_data))
model_other_unknown <- coeff_table(multi_lms_other_unknown)   %>%   
  mutate(clinical_var = "ART_other_unknown") %>% 
  select(community, clinical_var, everything())

# Adjust with FDR
all_stats_lipid_met_clinical <- rbind(model_Ethnicity, model_immunodecifiency, model_INSTI, model_PI, model_NNRTI, model_NRTI, model_other_unknown, model_SAT, model_thym_exposure, model_VAT) %>% 
  mutate(pvalue_adjusted = round(p.adjust(.$`Pr(>|t|)`, "fdr"),5))

# From character to numeric
cols = c("Estimate","Std. Error", "t value", "Pr(>|t|)", "pvalue_adjusted")    
all_stats_lipid_met_clinical[,cols] = apply(all_stats_lipid_met_clinical[,cols], 2, function(x) as.numeric(as.character(x)));

model2_all_stats_lipid_met_clinical <- cbind(all_stats_lipid_met_clinical[, c("community", "clinical_var", "Biochemicals")], 
             round(all_stats_lipid_met_clinical[,c("Estimate","Std. Error", "t value")],2), 
             round(all_stats_lipid_met_clinical[, c("Pr(>|t|)", "pvalue_adjusted")], 6)) %>% 
  write_csv("data/model2_lipid_met_vs_clinical.csv")
