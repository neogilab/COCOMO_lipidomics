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

# Change TAG names to accomodate the lipidomR function "map_lipid_names", insert paranthesis
##col_names_1 <- gsub("\\S[A-Z]*[0-9]*\\S[0-9]$", "", colnames(lipidomics_data_tidy))
##col_names_2 <- gsub('^(TAG)(.*)$', '\\1(\\2)', col_names_1)

col_names <- gsub('^(TAG)([0-9]*\\S[0-9]*)(\\SFA)([0-9]*\\S[0-9]*)$', '\\1(\\2)\\3(\\4)', colnames(lipidomics_data_tidy))
colnames(lipidomics_data_tidy) <- col_names

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

figure.output <-
  heatmap_lipidome_from_limma(
    x = result.limma$"model",
    names.mapping = names.mapping,
    axis.x.carbons = FALSE,
    class.facet = "row",
    plot.all = TRUE,
    plot.individual = FALSE,
    print.figure = TRUE,
    scales = "free",
    space = "free"
  )

# Create factor-specific figures.
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

# Print the figure of differences between ConditionHIV_NoMetS and ConditionHIV_MetS.
print( figure.output[[ "GroupHIV_MetS" ]] )

