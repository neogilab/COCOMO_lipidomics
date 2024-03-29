##########################################################
# Master script running all the analyses scripts
##########################################################

# Clear workspace
rm(list = ls())

# Run project
source("R/00_package_dependencies.R")
source("R/01_load_data.R")
source("R/02_tidy_data.R")
source("R/03_data_distribution.R")
source("R/04_unsupervised_analysis_lipidomics.R")
source("R/05_supervised_analysis_lipidomics.R") # OBS: Takes long time to run (> 1 hour, because of MUVR)
source("R/06_univariate_analysis_lipidomics.R")
source("R/07_consensus_analysis_lipidomics.R")
source("R/08_systematical_lipid_composition.R")
source("R/09_lipid_ontology_enrichment_analysis.R")
source("R/10_integration_clinical_and_omics_data.R")