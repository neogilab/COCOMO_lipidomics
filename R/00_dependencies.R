##########################################################
# Setup script with necessary packages
# Will be updated....
##########################################################

# Clear workspace
rm(list = ls())

# Used packages
list_of_packages <- c("readxl","tidyverse", "patchwork", "devtools", "ropls", "ggrepel", "limma")

# Install packages only if not previously installed
new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[, "Package"])]
if(length(new_packages)) {
  install.packages(new_packages)
}