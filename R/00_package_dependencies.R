##########################################################
# Setup script with necessary packages to run all scripts
##########################################################

# Clear workspace
rm(list = ls())

# Used packages
list_of_packages <- c("readxl","tidyverse", "caret", "patchwork", "Publish", "broom", "ropls", 
                      "randomForest", "devtools", "doParallel", "MUVR", "plotROC",
                      "limma", "UpSetR", "RVenn", "lipidomeR", "PCAtools")

# Install packages only if not previously installed
new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[, "Package"])]
if(length(new_packages)) {
  install.packages(new_packages)
}