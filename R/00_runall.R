##########################################################
# Master script running all the analyses scripts
##########################################################

# Clear workspace
rm(list = ls())

# Run project
source("R/00_dependencies.R")
source("R/01_load_data.R")
source("R/02_tidy_data.R")
source("R/03_data_distribution.R")
source("R/04_unsupervised_analysis.R")
source("R/05_supervised_analysis.R")
source("R/06_univariate_analysis.R")
source("R/07_consensus_analyses.R")
source("R/09_lipid_ontology_enrichment_analysis.R")
