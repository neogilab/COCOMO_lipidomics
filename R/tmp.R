##########################################################
# Temporary script, with experimental code from 
# 05_supervised_analysis.R
##########################################################

# Clear workspace
rm(list = ls())

library(readxl)
library(readr)
library(tidyverse)



# Random forest with randomForest -------------------------------------
# RF link: https://machinelearningmastery.com/tune-machine-learning-algorithms-in-r/
# Create features and target
X <- lipidomics_data_tidy_hiv_raw[, 3:ncol(lipidomics_data_tidy_hiv_raw)]
y <- as.factor(lipidomics_data_tidy_hiv_raw$Condition)

# Set seed for reproducibility
set.seed(seed)

################# Dividing data
# Split data randomly into training and validation sets
index <- createDataPartition(y, p=0.80, list=FALSE)
X_train <- X[ index, ]
X_val <- X[-index, ]
y_train <- y[index]
y_val <-y[-index]

################# Tuning the parameters and applying cross validation
# 10 folds repeat 3 times
control <- trainControl(method="repeatedcv", 
                        number=10, 
                        repeats=3)
metric <- "Accuracy" #Accuracy is the percentage of correctly classifies instances out of all instances. Useful for binary classification
method <- "rf"
set.seed(seed)

# Number randomely variable selected is mtry
mtry <- sqrt(ncol(X))
tunegrid <- expand.grid(.mtry=mtry)

# Training the model with train data
rf_gridsearch <- train(Condition~., 
                       data=lipidomics_data_tidy_hiv_raw[, 3:ncol(lipidomics_data_tidy_hiv)], 
                       method=method, 
                       metric=metric, 
                       tuneGrid=tunegrid, 
                       trControl=control,
                       importance=TRUE)

# Selecting the final/optimal rf model 
rf_model <- rf_gridsearch$finalModel

# Selecting the important/significant lipids
importance <- as.data.frame(rf_model$importance)
sign_lipids_HIV_randomForest <- subset(rownames(importance), importance$MeanDecreaseAccuracy >= acc_thres)
## rf_sig_lipids <- subset(rownames(importance), importance$MeanDecreaseGini >= gini_thres)


png("results/05_rf_feature_importance_acc_HIV.png")
varImpPlot(rf_model, main ='Feature importance', n.var = 15)
dev.off()

################# Making predictions with validation data and visualize
pred <- predict(rf_gridsearch, X_val)
result <- X_val
result['condition'] <- y_val
result['prediction'] <-  pred

# Performance of binary classification
cm <- confusionMatrix(table(y_val, pred), cutoff = 0.5)
cmtable <- cm$tabella

#Draw the ROC curve 
rf_probs <- predict(rf_gridsearch, X_val, type="prob")

# Plot The ROC curve
png("results/05_rf_ROCcurve_HIV.png")
rocplot <- ggplot(X_val, aes(m = rf_probs$HIV_NoMetS, d = y_val))+ geom_roc(n.cuts=20,labels=FALSE)
rocplot + 
  style_roc(theme = theme_grey) + 
  ggtitle("ROC plot for Random Forest model") +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = .75, y = .25, label = paste("AUC =", round(calc_auc(rocplot)["AUC"], 3)))
dev.off()


# Support Vector Machince (SVM) -------------------------------------------
# Create features and target
X <- lipidomics_data_tidy_hiv_raw[, 3:ncol(lipidomics_data_tidy_hiv_raw)]
y <- as.factor(lipidomics_data_tidy_hiv_raw$Condition)

# Split data randomly into training and validation sets
index <- createDataPartition(y, p=0.80, list=FALSE)
X_train <- X[ index, ]
X_val <- X[-index, ]
y_train <- y[index]
y_val <-y[-index]


# 10 folds repeat 3 times
control <- trainControl(method="repeatedcv", 
                        number=10, 
                        repeats=3,  
                        savePred=TRUE, 
                        classProb=TRUE)

# Fit the model on the training set
set.seed(seed)
svm_fit <- train(Condition~., 
                 data=lipidomics_data_tidy_hiv_raw[, 3:ncol(lipidomics_data_tidy_hiv)],
                 method = "svmRadial", #The radial basis kernel is extremely flexible and as a rule of thumb, we generally start with this kernel when fitting SVMs in practice.
                 trControl = control,
                 preProcess = c("center","scale"),
                 tuneLength = 10,
                 importance=TRUE
)


head(svm_fit$pred)

xx <- svm_fit$pred

# Selecting the final/optimal svm model 
svm_model <- svm_fit$finalModel



pred <- predict(svm_fit, X_val)
result <- X_val
result['condition'] <- y_val
result['prediction'] <-  pred



confusionMatrix(iris$Species, predict(m1))


# Make predictions on the test data
predicted.classes <- svm_fit %>% predict(X_val)
# Compute model accuracy rate
mean(predicted.classes == X_val$Condition)








# Orthogonal partial least-squares (OPLS-DA) ---------------------------------------------------------------------
# OPLS-DA on HIV groups
png("results/05_opls-da_HIV.png")
lipid_oplsda <- opls(lipidomics_data_tidy_hiv[, 4:ncol(lipidomics_data_tidy_hiv)], lipidomics_data_tidy_hiv$Condition,
                     predI = 1, orthoI = NA)
dev.off()
## When using subset the following error message: 
#No model was built because the first orthogonal component was already not significant;
#Select a number of orthogonal components of 1 if you want the algorithm to compute a model despite this.

# OPLS-DA on HIV groups: Extract VIP scores for all lipids, cut-off is usally >1 for vipVn values
lipid_oplsda_vipvn <- getVipVn(lipid_oplsda)
lipid_oplsda_vipvn <- as.data.frame(lipid_oplsda_vipvn)
lipid_oplsda_vipvn_df <- rownames_to_column(lipid_oplsda_vipvn, "Lipids")
lipid_oplsda_vipvn_df_order <- lipid_oplsda_vipvn_df[order(lipid_oplsda_vipvn, decreasing = TRUE),]

# OPLS-DA: Plot vip values
oplsda_vipplot <- lipid_oplsda_vipvn_df_order %>% 
  filter(!lipid_oplsda_vipvn <= vipvn_thres) %>% 
  ggplot(aes(x = reorder(Lipids, -lipid_oplsda_vipvn), y = lipid_oplsda_vipvn)) + 
  geom_col() + 
  labs(title = "OPLS-DA: Vip scores",
       subtitle = "Cut-off vip value < 1",
       x = "Lipids", 
       y = "Variable Importance in Projection ") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("results/05_vip_values_oplsda.png", plot = oplsda_vipplot, device = "png", width = 8, height = 3.1)

# OPLS-DA on HIV groups: Significant lipids
sign_lipids_HIV_oplsda <- subset(lipid_oplsda_vipvn_df_order$Lipids, lipid_oplsda_vipvn_df_order$lipid_oplsda_vipvn > vipvn_thres)



# biosigner (PLS-DA, Random Forest, SVM): HIVnomets vs. HIVmets --------------------
# Train and test
# Cite biosigner article: https://www.bioconductor.org/packages/devel/bioc/vignettes/biosigner/inst/doc/biosigner-vignette.html
trainVi <- 1:floor(0.8*nrow(lipidomics_data_tidy_hiv[, 4:ncol(lipidomics_data_tidy_hiv)]))
testVi <- setdiff(1:nrow(lipidomics_data_tidy_hiv[, 4:ncol(lipidomics_data_tidy_hiv)]), trainVi)

png("results/05_biosigner_plsda_rf_svm.png")
lipid_train_all <- biosign(x = lipidomics_data_tidy_hiv[trainVi, 4:ncol(lipidomics_data_tidy_hiv)], 
                           y = lipidomics_data_tidy_hiv[trainVi, c("Condition")],
                           bootI = 50, 
                           methodVc = "all")
dev.off()

# PLS-DA: Significant lipids
sign_lipids_HIV_biosign_plsda <- getSignatureLs(lipid_train_all)$plsda

# Random Forest: Significant lipids 
sign_lipids_HIV_biosign_rf <- getSignatureLs(lipid_train_all)$randomforest

# SVM: Significant lipids
sign_lipids_HIV_biosign_svm <- getSignatureLs(lipid_train_all)$svm


# fitted types
FitDF <- biosigner::predict(lipid_train_all)

# Confussion tables for each method
lapply(FitDF, function(predFc) table(actual = lipidomics_data_tidy_hiv[trainVi, c("Condition")], 
                                     predicted = predFc))

# Balanced accuracies
sapply(FitDF, function(predFc) {
  conf <- table(lipidomics_data_tidy_hiv[trainVi, c("Condition")], predFc)
  conf <- sweep(conf, 1, rowSums(conf), "/")
  round(mean(diag(conf)), 3)
})
# plsda         randomforest        svm 
# 0.828         0.757               0.789 


round(biosigner::getAccuracyMN(lipid_train_all)["S", ], 3)


# Performance of the test set
TestDF <- biosigner::predict(lipid_train_all, newdata = lipidomics_data_tidy_hiv[trainVi, 4:ncol(lipidomics_data_tidy_hiv)])
sapply(TestDF, function(predFc) {
  conf <- table(lipidomics_data_tidy_hiv[trainVi, c("Condition")], predFc)
  conf <- sweep(conf, 1, rowSums(conf), "/")
  round(mean(diag(conf)), 3)
})
#  plsda        randomforest          svm 
#  0.831        1.000                 0.818 


# List with significant lipids from different methods -------------------------------
#method_list_supervised <- list('ropls: PLS-DA' = sign_lipids_HIV_plsda,
# 'ropls: OPLS-DA' = sign_lipids_HIV_oplsda,
#  'Biosign: PLS-DA' = sign_lipids_HIV_biosign_plsda,
# 'Biosign: Random Forests' = sign_lipids_HIV_biosign_rf,
#  'Biosign: SVM' = sign_lipids_HIV_biosign_svm,
# 'MUVR: Random Forest' = sign_lipids_HIV_MUVR)

#method_list_supervised <- list('PLS-DA with "ropls"' = sign_lipids_HIV_plsda,
#                              'Random Forest with "MUVR"' = sign_lipids_HIV_MUVR,
#                             'Random Forest with "randomForest"' = sign_lipids_HIV_randomForest)

method_list_supervised <- list('PLS-DA with "ropls"' = sign_lipids_HIV_plsda,
                               'Random Forest with "MUVR"' = sign_lipids_HIV_MUVR)

# Save list to a file
save(method_list_supervised, file="data/05_methods_supervised.RData")