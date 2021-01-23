##########################################################
# Script identifying the trends in the systamtical 
# Composition of the lipids from lipidomics data 
##########################################################

# Clear workspace
rm(list = ls())

# Load packages
library(lipidomeR)
library(tidyverse)

# Load data
lipidomics_data_tidy <- read_csv("data/02_lipidomics_data_tidy.csv") %>% 
  mutate(Group = case_when(Condition == "HIV_NoMetS" ~ "HIV_Control", 
                           Condition == "HIV_MetS" ~ "HIV_MetS")) %>% 
  select(Group, everything(), -c("GENDER","ID_Condition", "Condition")) 



# Limma model -----------------------------------------------------------------
# List of lipid names
lip_names <- colnames(lipidomics_data_tidy[,2:ncol(lipidomics_data_tidy)])

# Create a mapping of the lipid names.
names.mapping <- map_lipid_names(x = lip_names)

# Compute the regression models
result.limma <-
  compute_models_with_limma(
    x = lipidomics_data_tidy, #[1:400],
    dependent.variables = names.mapping$"Name",
    independent.variables = c("Group"))



# Heatmap of lipid classes -----------------------------------------------------------------
# Create factor-specific figures
figure.output <-
  heatmap_lipidome_from_limma(
    x = result.limma$"model",
    names.mapping = names.mapping,
    axis.x.carbons = FALSE,
    class.facet = "wrap",
    omit.class = "PA",
    plot.all = FALSE,
    plot.individual = TRUE,
    print.figure = FALSE,
    scales = "free",
    space = "free"
  )

# Save the figure of differences between Condition[HIV_NoMet/HIV_MetS]
ggsave("results/08_heatmaps_lipidomeR.png", plot = figure.output[[ "GroupHIV_MetS" ]], device = "png", height = 11) #, width = 6.17,