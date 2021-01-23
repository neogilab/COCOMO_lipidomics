##########################################################
# Script for supervised analysis of lipidomics data
# This script runs dimensionality reduction by PLS-DA and
# Classification models such as logistic regression and
# Random Forests. Performance are evaluted by Q2 and ROC/AUC
##########################################################

# Clear workspace
rm(list = ls())

# Load packages
library(tidyverse)
library(ropls)
library(randomForest)
library(devtools)       
library(doParallel)     # Parallel processing (install.packages("doParallel", repos="http://R-Forge.R-project.org") )
library(MUVR)           # Multivariate modelling (install_git("https://gitlab.com/CarlBrunius/MUVR.git"))
library(plotROC)    

# Load data 
lipidomics_data_tidy <- read_csv("data/02_lipidomics_data_tidy.csv")  

# Set a seed for reproducible data. Seed sampled with -> sample(1e6, 1)
set.seed(878231)

# Variables
vipvn_thres <- 1 # Cut-off is usally > 1 for vipVn values



# Log transform data -----------------------------------------------------------
lipidomics_data_tidy[, 4:ncol(lipidomics_data_tidy)] <- sapply(log2(lipidomics_data_tidy[, 4:ncol(lipidomics_data_tidy)]), as.numeric)
lipidomics_data_tidy <- as.data.frame(lipidomics_data_tidy)
rownames(lipidomics_data_tidy) <- lipidomics_data_tidy$ID_Condition



# Partial least-squares (PLS-DA) ---------------------------------------------------------------------
# PLS-DA on Condition[HIV-NoMetS/HIV-MetS]
png("results/05_1_plsda_condition.png")
plsda_condition <- opls(lipidomics_data_tidy[, 4:ncol(lipidomics_data_tidy)], lipidomics_data_tidy$Condition)
dev.off()

# PLS-DA on Sex[Male/Female]
png("results/05_2_plsda_sex.png")
plsda_gender <- opls(lipidomics_data_tidy[, 4:ncol(lipidomics_data_tidy)], lipidomics_data_tidy$GENDER)
dev.off()

# Extract VIP scores (PLS-DA on Condition[HIV-NoMetS/HIV-MetS])
plsda_condition_vipvn <- getVipVn(plsda_condition)
plsda_condition_vipvn <- as.data.frame(plsda_condition_vipvn)
plsda_condition_vipvn_df <- rownames_to_column(plsda_condition_vipvn, "Lipids")
plsda_condition_vipvn_df_order <- plsda_condition_vipvn_df[order(plsda_condition_vipvn, decreasing = TRUE),]

# Significant lipids (PLS-DA on Condition[HIV-NoMetS/HIV-MetS])
sign_lipids_HIV_plsda <- subset(plsda_condition_vipvn_df_order$Lipids, plsda_condition_vipvn_df_order$plsda_condition_vipvn > vipvn_thres)



# Random forest with MVUR -------------------------------------
# MVUR = Multivariate methods with Unbiased Variable selection -> Method used when few observations and a large number of variables

# Create "condition" factor with two levels "HIV_MetS" and "HIV_NoMetS"
condition_fac <- as.factor(lipidomics_data_tidy$Condition)

# Set method parameters, for parallel processing
nCore=detectCores()-1   # Number of processor threads to use, uses all but one thread
nRep=40              # Number of MUVR repetitions, usally between 20 and 50 
nOuter=8                # Number of outer cross-validation segments, usally between 6 and 8. Higher number when fewer observations ???> increase number of observations in the model training
varRatio=0.85           # Proportion of variables kept per iteration, usally start out low 0.75 and increase towards 0.85-0.9 for final processing 
method='RF'             # Selected core modelling algorithm to Random Forest

# OBS: hard to set seed for reproducibility, as it is parallel processing 

# Set up parallel processing using doParallel 
cl=makeCluster(nCore)   
registerDoParallel(cl)

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
  geom_point(stat='identity', aes(col=`Lipid class`), size=5)  +
  scale_color_manual(values = c("#32CD32", "#6495ED")) +
  labs(title="VIP score plot") + 
  ylab("Rank") + xlab("Increasing importance to group separation") +
  ylim(-50, 500) +
  coord_flip() +
  geom_segment(aes(x=0, xend = 13.5 , y=-50, yend = -50), size=1, arrow = arrow(length = unit(0.3,"cm")))
ggsave("results/05_3_vip_score_plot_MUVR.png", plot = vip_score_plot, device = "png", width = 6.17, height = 5)



# Performance and accuracy of models --------------------------------------
#  PLS-DA -> Evaluate Q2 value
plsda_condition

# ROC curve MUVR: "Min" + "Mid" + "Max" model
png("results/05_4_roc_MUVR_MinMidMax.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 13)        # font size)     
par(pty="s")
roc(classModel$inData$Y, 
    classModel$yPred$min[,1], 
    plot = TRUE, 
    legacy.axes = TRUE,
    percent = TRUE,
    xlab = "False Positive Percentage",
    ylab = "True Positive Percentage",
    col = "#4daf4a", #"#377eb8",
    lwd = 3,
    print.auc = TRUE,
    print.auc.x = 40,
    print.auc.y = 50,
    levels=c("HIV_NoMetS", "HIV_MetS")) 
plot.roc(classModel$inData$Y, 
         classModel$yPred$mid[,1], 
         percent = TRUE,
         col = "#377eb8",
         lwd = 3,
         print.auc = TRUE,
         print.auc.x = 40,
         print.auc.y = 45,
         add = TRUE,
         levels=c("HIV_NoMetS", "HIV_MetS")) 
plot.roc(classModel$inData$Y, 
         classModel$yPred$max[,1], 
         percent = TRUE,
         col = "red",
         lwd = 3,
         print.auc = TRUE,
         print.auc.x = 40,
         print.auc.y = 40,
         add = TRUE,
         levels=c("HIV_NoMetS", "HIV_MetS")) 
legend("bottomright", legend = c("'Min' MUVR model", "'Mid' MUVR model", "'Max' MUVR model"), col = c("#4daf4a", "#377eb8", "red"), lwd = 3)
par(pty="m")
dev.off()



# List with significant lipids from different methods -------------------------------
method_list_supervised <- list('PLS-DA with "ropls"' = sign_lipids_HIV_plsda,
                              'Random Forest with "MUVR"' = sign_lipids_HIV_MUVR)

# Save list to a file
save(method_list_supervised, file="data/05_methods_supervised.RData")