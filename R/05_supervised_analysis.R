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
library(plotROC)    

# Variables
vipvn_thres <- 1 # Cut-off is usally > 1 for vipVn values


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
##lipidomics_data_tidy_hiv <- lipidomics_data_tidy %>% 
  ##filter(!str_detect(GENDER, "Female")) 



# Partial least-squares (PLS-DA) ---------------------------------------------------------------------
# PLS-DA on three groups
png("results/05_pls-da_condition.png")
plsda_condition <- opls(lipidomics_data_tidy[, 4:ncol(lipidomics_data_tidy)], lipidomics_data_tidy$Condition)
dev.off()

# PLS-DA on HIV patients divided by condition
png("results/05_pls-da_gender.png")
plsda_gender <- opls(lipidomics_data_tidy[, 4:ncol(lipidomics_data_tidy)], lipidomics_data_tidy$GENDER)
dev.off()

# PLS-DA on HIV patients divided by gender
#png("results/05_pls-da_HIV_gender.png")
#lipid_plsda <- opls(lipidomics_data_tidy_hiv[, 4:ncol(lipidomics_data_tidy_hiv)], lipidomics_data_tidy_hiv$GENDER)
#dev.off()

# PLS-DA on HIV groups: Extract VIP scores for all lipids
plsda_condition_vipvn <- getVipVn(plsda_condition)
plsda_condition_vipvn <- as.data.frame(plsda_condition_vipvn)
plsda_condition_vipvn_df <- rownames_to_column(plsda_condition_vipvn, "Lipids")
plsda_condition_vipvn_df_order <- plsda_condition_vipvn_df[order(plsda_condition_vipvn, decreasing = TRUE),]

# PLS-DA: Plot vip values
plsda_vipplot <- plsda_condition_vipvn_df_order %>% 
  filter(!plsda_condition_vipvn <= vipvn_thres) %>% 
  ggplot(aes(x = reorder(Lipids, -plsda_condition_vipvn), y = plsda_condition_vipvn)) + 
  geom_col() + 
  labs(title = "PLS-DA: Vip scores",
       subtitle = "Cut-off vip value < 1",
       x = "Lipids", 
       y = "Variable Importance in Projection ") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("results/05_vip_values_plsda.png", plot = plsda_vipplot, device = "png", width = 8, height = 3.1)

# PLS-DA on HIV groups: Significant lipids
sign_lipids_HIV_plsda <- subset(plsda_condition_vipvn_df_order$Lipids, plsda_condition_vipvn_df_order$plsda_condition_vipvn > vipvn_thres)


# Random forest with MVUR -------------------------------------
# MVUR = Multivariate methods with Unbiased Variable selection
# Method used when few observations and a large number of variables

# Filter out the control group
##lipidomics_data_tidy_HIV <- lipidomics_data_tidy %>% 
##  filter(!str_detect(Condition, "Ctrl")) 

# Create "condition" factor with two levels "HIV_MetS" and "HIV_NoMetS"
condition_fac <- as.factor(lipidomics_data_tidy$Condition)

# Set method parameters, for parallel processing
nCore=detectCores()-1   # Number of processor threads to use, uses all but one thread
nRep=50              # Number of MUVR repetitions, usally between 20 and 50 
nOuter=8                # Number of outer cross-validation segments, usally between 6 and 8. Higher number when fewer observations ???> increase number of observations in the model training
varRatio=0.85           # Proportion of variables kept per iteration, usally start out low 0.75 and increase towards 0.85-0.9 for final processing 
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
plotMV(classModel, model='max')         # Look at the model of choice: min, mid or max
plotStability(classModel, model='max')  # The stability plot for classification analysis generates three subplots. 1. Number of selected variables for each repetition as well as cumulative average over the repetitions; 2. The proportion of selected variables reports the ratio of the final variable selection found in each repetition and cumulatively, averaged over the number of repetitions; 3. Number of misclassifications per repetition and cumulatively. 
plotVIP(classModel, model='max')        # Boxplot of the variables automatically selected from optimal modelling performance. 
getVIP(classModel, model='max')         # Extract most informative variables: Lower rank is better


# Extract significant lipids
rf_MUVR_sig <- getVIP(classModel, model='max') %>% 
  rename(lipids = name) %>% 
  mutate(`Lipid class` = case_when(str_detect(lipids, "DAG") ~ "Diacylglycerol",
                                   str_detect(lipids, "TAG") ~ "Triacylglycerol"))
sign_lipids_HIV_MUVR <- rf_MUVR_sig$lipids

# order
rf_MUVR_sig$lipids <- factor(rf_MUVR_sig$lipids, levels = rf_MUVR_sig$lipids[order(desc(rf_MUVR_sig$rank))])

# Plot significant lipids according to their rank
vip_score_plot <- ggplot(as.data.frame(rf_MUVR_sig), aes(x=lipids, y=rank, label=round(rank,1))) + 
  geom_point(stat='identity', aes(col=`Lipid class`), size=4)  +
 # scale_color_manual(labels = c("Triacylglycerol", "Diacylglycerol"), 
  #                   values = c("Triacylglycerol"="#00ba38", "Diacylglycerol"="#f8766d")) + 
  geom_text(color="black", size=2) +
  labs(title="VIP score plot", 
       subtitle="Lower rank indicates better prediction variables") + 
  ylab("Rank") + xlab("Increasing importance to group seperation") +
  ylim(-50, 470) +
  coord_flip() +
  geom_segment(aes(x=0, xend = 13.5 , y=-50, yend = -50), size=1, arrow = arrow(length = unit(0.3,"cm")))
vip_score_plot

ggsave("results/05_vip_score_plot.png", plot = vip_score_plot, device = "png")#, width = 6.17, height = 3.1)


# Performance of binary classification
cm <- classModel$Fit$rfFitMid$confusion

#Draw the ROC curve ..
MUVR_probs <- as.data.frame(classModel$yPred$max)
auc_mid_MUVR <- classModel$auc[3,1]

roc_data <- cbind(classModel$yClass$max, MUVR_probs) %>% 
  mutate(binary_condition = case_when(classModel$yClass$mid == 'HIV_NoMetS' ~ 0,
                                      classModel$yClass$mid == 'HIV_MetS' ~ 1))

rocplot <- ggplot(lipidomics_data_tidy[, 4:ncol(lipidomics_data_tidy)], 
                  aes(m = roc_data$HIV_MetS, 
                      d = roc_data$binary_condition)) + 
  geom_roc(n.cuts=20) + 
  style_roc(theme = theme_grey) + 
  ggtitle("ROC plot for Random Forest model MVUR") +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = .75, y = .25, label = paste("AUC =", round(auc_mid_MUVR, 3)))
rocplot
ggsave("results/05_roc.png", plot = rocplot, device = "png")#, width = 6.17, height = 3.1)



# List with significant lipids from different methods -------------------------------
method_list_supervised <- list('PLS-DA with "ropls"' = sign_lipids_HIV_plsda,
                              'Random Forest with "MUVR"' = sign_lipids_HIV_MUVR)

# Save list to a file
save(method_list_supervised, file="data/05_methods_supervised.RData")

