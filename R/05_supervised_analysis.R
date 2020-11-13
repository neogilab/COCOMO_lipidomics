##########################################################
# Script for supervised analysis of lipidomics data
# This script runs dimensionality reduction by PLS and PLS-DA
# Furthermore classification models such as logistic regression
# and Random Forests are performed 
##########################################################

# Clear workspace
rm(list = ls())

# Load packages
library(tidyverse)
library(ropls)
library(biosigner)
library(randomForest)
library(devtools)       # For MVUR
library(doParallel)     # Parallel processing (install.packages("doParallel", repos="http://R-Forge.R-project.org") )
library(MUVR)           # Multivariate modelling (install_git("https://gitlab.com/CarlBrunius/MUVR.git"))
library(caret)          # For function "createDataPartition"
library(plotROC)    

# Variables
vipvn_thres <- 1 # Cut-off is usally > 1 for vipVn values
acc_thres <- 0.0025 
gini_thres <- 0.5



# Load data ---------------------------------------------------------------
lipidomics_data_tidy <- read_csv("data/02_lipidomics_data_tidy.csv")  

# Filter out the control group
#lipidomics_data_tidy_hiv_raw <- lipidomics_data_tidy %>% 
  #filter(!str_detect(GENDER, "Female")) 



# Set seed ----------------------------------------------------------------
# Set a seed in order to be able to create reproducible data. Sample seed by the function: sample(1e6, 1)
seed <- 878231



# Transform data -----------------------------------------------------------
# Log transform data
lipidomics_data_tidy[, 4:ncol(lipidomics_data_tidy)] <- sapply(log2(lipidomics_data_tidy[, 4:ncol(lipidomics_data_tidy)]), as.numeric)

# Transform data
lipidomics_data_tidy <- as.data.frame(lipidomics_data_tidy)
rownames(lipidomics_data_tidy) <- lipidomics_data_tidy$ID_Condition

# Filter out the females among HIV infected
lipidomics_data_tidy_hiv <- lipidomics_data_tidy %>% 
  filter(!str_detect(GENDER, "Female")) 



# Partial least-squares (PLS-DA) ---------------------------------------------------------------------
# PLS-DA on three groups
png("results/05_pls-da.png")
opls(lipidomics_data_tidy[, 4:ncol(lipidomics_data_tidy_hiv)], lipidomics_data_tidy$Condition)
dev.off()

# PLS-DA on HIV patients divided by condition
png("results/05_pls-da_HIV_condition.png")
lipid_plsda <- opls(lipidomics_data_tidy_hiv[, 4:ncol(lipidomics_data_tidy_hiv)], lipidomics_data_tidy_hiv$Condition)
dev.off()

# PLS-DA on HIV patients divided by gender
#png("results/05_pls-da_HIV_gender.png")
#lipid_plsda <- opls(lipidomics_data_tidy_hiv[, 4:ncol(lipidomics_data_tidy_hiv)], lipidomics_data_tidy_hiv$GENDER)
#dev.off()

# PLS-DA on HIV groups: Extract VIP scores for all lipids
lipid_plsda_vipvn <- getVipVn(lipid_plsda)
lipid_plsda_vipvn <- as.data.frame(lipid_plsda_vipvn)
lipid_plsda_vipvn_df <- rownames_to_column(lipid_plsda_vipvn, "Lipids")
lipid_plsda_vipvn_df_order <- lipid_plsda_vipvn_df[order(lipid_plsda_vipvn, decreasing = TRUE),]

# PLS-DA: Plot vip values
plsda_vipplot <- lipid_plsda_vipvn_df_order %>% 
  filter(!lipid_plsda_vipvn <= vipvn_thres) %>% 
  ggplot(aes(x = reorder(Lipids, -lipid_plsda_vipvn), y = lipid_plsda_vipvn)) + 
  geom_col() + 
  labs(title = "PLS-DA: Vip scores",
       subtitle = "Cut-off vip value < 1",
       x = "Lipids", 
       y = "Variable Importance in Projection ") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("results/05_vip_values_plsda.png", plot = plsda_vipplot, device = "png", width = 8, height = 3.1)

# PLS-DA on HIV groups: Significant lipids
sign_lipids_HIV_plsda <- subset(lipid_plsda_vipvn_df_order$Lipids, lipid_plsda_vipvn_df_order$lipid_plsda_vipvn > vipvn_thres)


# Random forest with MVUR -------------------------------------
# MVUR = Multivariate methods with Unbiased Variable selection
# Method used when few observations and a large number of variables

# Create "condition" factor with two levels "HIV_MetS" and "HIV_NoMetS"
condition_fac <- as.factor(lipidomics_data_tidy$Condition)

# Set method parameters, for parallel processing
nCore=detectCores()-1   # Number of processor threads to use, uses all but one thread
nRep=nCore              # Number of MUVR repetitions, usally between 20 and 50 
nOuter=8                # Number of outer cross-validation segments, usally between 6 and 8. Higher number when fewer observations –> increase number of observations in the model training
varRatio=0.8            # Proportion of variables kept per iteration, usally start out low 0.75 and increase towards 0.85-0.9 for final processing 
method='RF'             # Selected core modelling algorithm to Random Forest

# Set up parallel processing using doParallel 
cl=makeCluster(nCore)   
registerDoParallel(cl)

# OBS: hard to set seed for reproducibility, as it is parallel processing 

# Perform modelling
classModel = MUVR(X=lipidomics_data_tidy[, 4:ncol(lipidomics_data_tidy)], 
                  Y=condition_fac, 
                  nRep=nRep, 
                  nOuter=nOuter, 
                  varRatio=varRatio, 
                  method=method,
                  importance=TRUE)

# Stop parallel processing
stopCluster(cl)
  
# Examine model performance and output
classModel$miss                         # Number of misclassifications for min, mid and max models
classModel$nVar                         # Number of variables for min, mid and max models
cbind(condition_fac, classModel$yClass)          # Actual class side-by-side with min, mid and max predictions
plotVAL(classModel)
plotMV(classModel, model='min')         # Look at the model of choice: min, mid or max
plotStability(classModel, model='min')  # The stability plot for classification analysis generates three subplots. 1. Number of selected variables for each repetition as well as cumulative average over the repetitions; 2. The proportion of selected variables reports the ratio of the final variable selection found in each repetition and cumulatively, averaged over the number of repetitions; 3. Number of misclassifications per repetition and cumulatively. 
plotVIP(classModel, model='min')        # Boxplot of the variables automatically selected from optimal modelling performance. 
getVIP(classModel, model='mid')         # Extract most informative variables: Lower rank is better

# Extract significant lipids
rf_MUVR_sig <- getVIP(classModel, model='mid')
sign_lipids_HIV_MUVR <- rf_MUVR_sig$name

# Performance of binary classification
cm <- classModel$Fit$rfFitMid$confusion

#Draw the ROC curve ..
##MUVR_probs <- as.data.frame(classModel$yPred$mid)
##auc_mid_MUVR <- classModel$auc[2,1]
##rocplot <- ggplot(lipidomics_data_tidy_hiv_raw[, 4:ncol(lipidomics_data_tidy_hiv_raw)], 
  #                aes(m = MUVR_probs$HIV_NoMetS, 
  #                    d = classModel$yClass$mid)) + 
  #geom_roc(n.cuts=20) + 
  #style_roc(theme = theme_grey) + 
  #ggtitle("ROC plot for Random Forest model MVUR") +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #annotate("text", x = .75, y = .25, label = paste("AUC =", round(classModel$auc[2,1], 3)))
##rocplot



# Linear regression -----------------------------------------------------
xxx


# List with significant lipids from different methods -------------------------------
method_list_supervised <- list('PLS-DA with "ropls"' = sign_lipids_HIV_plsda,
                              'Random Forest with "MUVR"' = sign_lipids_HIV_MUVR)

# Save list to a file
save(method_list_supervised, file="data/05_methods_supervised.RData")

