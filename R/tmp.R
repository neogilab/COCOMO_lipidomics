# Clear workspace
rm(list = ls())

# Load packages
library(lipidomeR)
library(tidyverse)

# Load data
network_data <- read_csv("data/08_cytoscape_pos_table_FDR.csv")


# Change TAG names to accomodate the lipidomR function "map_lipid_names"
bio_names_1 <- gsub("\\S[A-Z]*[0-9]*\\S[0-9]$", "", network_data$Biochemicals)
bio_names_2 <- gsub('^(TAG)(.*)$', '\\1(\\2)', bio_names_1)

network_data <- network_data %>% 
  mutate(Bionames = bio_names_2)  

  
c1_names <- network_data %>% 
  filter(community == "c1") %>% 
  select(Bionames) %>% 
  write_csv("c1_new_names.csv")

c2_names <- network_data %>% 
  filter(community == "c2") %>% 
  select(Bionames) %>% 
  write_csv("c2_new_names.csv")

c3_names <- network_data %>% 
  filter(community == "c3") %>% 
  select(Bionames) %>% 
  write_csv("c3_new_names.csv")

all_names <- network_data %>% 
  select(Bionames) %>% 
  write_csv("all_new_names.csv")
