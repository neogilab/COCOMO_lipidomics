##########################################################
# Script loading the raw data files, saving them to csv
##########################################################

# Clear workspace
rm(list = ls())

# Load packages
library(readxl)
library(tidyverse)


# Load clinical data from excel file --------------------------------------
# This dataset is a combination of the COCOMO_HIV_data.xlsx and the COCOMO_validated_data.xlsx
clinical_data <- read_xlsx("data/_raw/clinical_data_COCOMO.xlsx") %>% 
  write_csv("data/01_clinical_data.csv")

# Data from the Visceral Adipose Tissue (VAT) and Subcutaneous Abdominal Tissue (SAT)
vatsat_data <- read_xlsx("data/_raw/Fat_data.xlsx") %>% 
  write_csv("data/01_vatsat_data.csv")


# Load lipidomics data from excel file ----------------------------------------------------
lipidomics_data <- read_xlsx("data/_raw/clp_lipidomics_data_clean.xlsx") %>% 
  write_csv("data/01_lipidomics_data.csv")


# Load metabolomics data from excel file ----------------------------------------------------
metabolomics_data <- read_xlsx("data/_raw/hd4_metabolomics_data.xlsx", sheet = 3) %>%
  write_csv("data/01_metabolomics_data.csv")

