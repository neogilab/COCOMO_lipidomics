##########################################################
# Script for cleaning data. This script will tidy the data 
# and handle missing values and low variance.
##########################################################

# Clear workspace
rm(list = ls())

# Load packages
library(tidyverse)
library(caret)

# Load data ---------------------------------------------------------------
clinical_data <- read_csv("data/01_clinical_data.csv")
lipidomics_data <- read_csv("data/01_lipidomics_data.csv")
metabolomics_data <- read_csv("data/01_metabolomics_data.csv")



# Tidying clinical data  --------------------------------------------------
# More to come...
# Combine the "ID" col with the "Condition" col
clinical_data_modified <- as_tibble(clinical_data) %>%
  select(SUBJECT_ID, CONDITION, everything()) %>% 
  unite(col = "ID_Condition", "SUBJECT_ID":"CONDITION", sep = "_", remove = FALSE, na.rm = TRUE)



# Tidying lipidomics data ---------------------------------------------
# Remove rows as they contain clinical data, except subject_id
lipidomics_data_modified <- lipidomics_data %>% 
  slice(13:n()) 

# Save list with lipid info
lip_info <- lipidomics_data_modified %>% 
  select(BIOCHEMICAL, SUPER_PATHWAY,SUB_PATHWAY,	KEGG,	`Group   HMDB_ID`) %>% 
  slice(3:n()) %>% 
  rename(Biochemical = BIOCHEMICAL, Super_pathway = SUPER_PATHWAY, Sub_pathway = SUB_PATHWAY, KEGG_ID = KEGG, HMDB_ID = `Group   HMDB_ID`) %>% 
  filter(!str_detect(Biochemical, 'Total'))  %>% 
  mutate(Group = "lipid")

# Transpose lipid matrix and combine "ID" col with "Condition" col
lipidomics_data_modified <- lipidomics_data_modified %>% 
  rownames_to_column(var = "rowname") %>%
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from = rowname, values_from = value) %>% 
  select(!name) %>% 
  unite(col = "ID_Condition", "1":"2", sep = "_", remove = FALSE, na.rm = TRUE)

# Add lipid names to the header and rename the first columns
colnames(lipidomics_data_modified) <- slice(lipidomics_data_modified, 2)
colnames(lipidomics_data_modified)[1:3] <- c("ID_Condition", "ID", "Condition")

# Remove rows containing additional information on lipids & columns with total lipid class concentration
lipidomics_data_modified <- lipidomics_data_modified %>% 
  slice(10:n()) %>% 
  select(everything(), -starts_with("Total")) 

# Select male samples only and deselect ctrl group
lipidomics_data_modified <- merge(x = lipidomics_data_modified, y = clinical_data_modified[ , c("ID_Condition", "GENDER")], by = "ID_Condition", all.x=TRUE) %>% 
#  filter(!GENDER == "Female") %>% 
  filter(!Condition == "Ctrl") %>% 
  select(GENDER, everything(), -ID)

# Make the chemical columns numeric
lipidomics_data_modified[, 4:ncol(lipidomics_data_modified)] <- sapply(lipidomics_data_modified[, 4:ncol(lipidomics_data_modified)], as.numeric)

# Check for missing values in lipid matrix
any(is.na(lipidomics_data_modified)) #FALSE = No missing values

# Remove columns with zero variance or near zero variance
compute_nzv <- nearZeroVar(lipidomics_data_modified[, 4:ncol(lipidomics_data_modified)], 
                           saveMetrics = TRUE)

nzv_lipids <- compute_nzv %>% 
  filter(nzv == FALSE) %>% 
  rownames(.)

lipidomics_data_modified_var <- lipidomics_data_modified %>% 
  select(GENDER, ID_Condition, Condition, all_of(nzv_lipids))


# Tidying metabolomics data ---------------------------------------------
# Remove rows as they contain clinical data, except subject_id
metabolomics_data_modified <- metabolomics_data %>% 
  slice(16:n()) 

# List with significant metabolites found in prev. study [ref: metabolomics paper], 11 in total
sign_met <- c('1-carboxyethylleucine', '4-cholesten-3-one', '4-hydroxyglutamate', 'alpha-ketoglutarate', 'carotene diol (2)', 
              'gamma-glutamylglutamate', 'glutamate', 'glycerate', 'isoleucine', 'palmitoyl-sphingosine-phosphoethanolamine (d18:1/16:0)',
              'pimeloylcarnitine/3-methyladipoylcarnitine (C7-DC)')

# Save list with metabolite info
met_info <- metabolomics_data_modified %>% 
  slice(2:n()) %>% 
  select(...2, ...12, `CLIENT IDENTIFIER`, ...3, ...4) %>% # Biochemical, Kegg entry and HMDB entry 
  rename(Biochemical = ...2, KEGG_ID = ...12, HMDB_ID = `CLIENT IDENTIFIER`, Super_pathway = ...3, Sub_pathway = ...4) %>% 
  filter(Biochemical %in% sign_met)  %>% 
  mutate(Group = "metabolite")

# Transpose metabolite matrix and combine the "ID" col with the "Condition" col
metabolomics_data_modified <- metabolomics_data_modified %>% 
  select(everything(), -c(...1, ...3:`CLIENT IDENTIFIER`)) %>% 
  rownames_to_column(var = "rowname") %>%
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from = rowname, values_from = value) %>% 
  unite(col = "ID_Condition", "name":"1", sep = "_", remove = FALSE, na.rm = TRUE) 

# Add metabolite names to the header and rename the first columns
colnames(metabolomics_data_modified) <- slice(metabolomics_data_modified, 1)
colnames(metabolomics_data_modified)[1:3] <- c("ID_Condition", "ID", "Condition")

# Remove rows containing additional information on the chemicals & columns with the total lipid class concentration
metabolomics_data_modified <- metabolomics_data_modified %>% 
  slice(2:n()) 

# Select male samples only
# metabolomics_data_modified <- merge(x = metabolomics_data_modified, y = clinical_data_modified[ , c("ID_Condition", "GENDER")], by = "ID_Condition", all.x=TRUE) %>% 
#filter(!GENDER == "Female") %>% 
#  filter(!Condition == "Ctrl") %>% 
#  select(GENDER, everything(), -ID)

# Make the metabolite columns numeric
metabolomics_data_modified[, 4:ncol(metabolomics_data_modified)] <- sapply(metabolomics_data_modified[, 4:ncol(metabolomics_data_modified)], as.numeric)

# Check for missing values in the metabolite matrix
any(is.na(metabolomics_data_modified)) #FALSE = No missing values

# Remove columns with zero variance or near zero variance
compute_nzv <- nearZeroVar(metabolomics_data_modified[, 4:ncol(metabolomics_data_modified)], saveMetrics = TRUE)

nzv_metabolites <- compute_nzv %>% 
  filter(nzv == FALSE) %>% 
  rownames(.)

metabolomics_data_modified_var <- metabolomics_data_modified %>% 
  select(GENDER, ID_Condition, Condition, all_of(nzv_metabolites))

# Select significant metabolites
metabolomics_sign_met <- metabolomics_data_modified %>% 
  select(ID_Condition, sign_met)



# Merge lipidomics and metabolomics data ----------------------------------
lip_met_merge_data_entries <- rbind(met_info, lip_info)

# Note: Make function out of this
# Check for duplicates in HMDB
n_occur_hmdb <- data.frame(table(lip_met_merge_data_entries$HMDB_ID))
n_occur_hmdb[n_occur_hmdb$Freq > 1,]
dupl_hmdb <- lip_met_merge_data_entries[lip_met_merge_data_entries$HMDB_ID %in% n_occur_hmdb$Var1[n_occur_hmdb$Freq > 1],]

# Remove one of the duplicate entries
dat <- dupl_hmdb$HMDB_ID
toDelete <- seq(1, length(dat), 2)
toDelete <- dupl_hmdb[-toDelete, ]
lip_met_merge_data_entries <- lip_met_merge_data_entries %>% 
  filter(!Biochemicals %in% toDelete$Biochemicals)

# Remove biochemical if no KEGG or HMDB entry
lip_met_merge_data_entries <- lip_met_merge_data_entries[rowSums(is.na(lip_met_merge_data_entries[,2:3])) != ncol(lip_met_merge_data_entries[,2:3]), ]

# Merge the two data frames
lip_met_merge_data_var <- merge(x = lipidomics_data_modified_var, 
                            y = metabolomics_data_modified_var,
                            by = "ID_Condition",
                            all.x = TRUE) 
  
# Data which have a KEGG entry, HMDB entry or both
lip_met_merge_data_clean <- lip_met_merge_data_var[colnames(lip_met_merge_data_var) %in% lip_met_merge_data_entries$Biochemical] %>% 
  mutate(ID_Condition = lip_met_merge_data$ID_Condition, 
         Condition = lip_met_merge_data$Condition.x, 
         Gender = lip_met_merge_data$GENDER) %>% 
  select(ID_Condition, Condition, Gender, everything())

  

# Save data to csv  -------------------------------------------
lipidomics_data_modified_var %>% 
  write_csv("data/02_lipidomics_data_tidy.csv")

lip_info %>% 
  filter(Biochemical %in% nzv_lipids) %>% 
  write_csv("data/02_lipidomics_data_info.csv")

metabolomics_data_modified_var %>% 
  write_csv("data/02_metabolomics_data_tidy.csv")

metabolomics_sign_met %>% 
  write_csv("data/02_metabolomics_sign_metabolites.csv")

met_info %>% 
  write_csv("data/02_metabolomics_data_info.csv")

lip_met_merge_data_clean %>% 
  write_csv("data/02_metabolomics_lipidomics_data_tidy.csv")

lip_met_merge_data_entries %>% 
  write_csv("data/02_metabolomics_lipidomics_data_info.csv")

clinical_data_modified %>% 
  write_csv("data/02_clinical_data_tidy.csv")







# Lipid class concentration -----------------------------------------------
### slet evt.
xxxx <- aggregate(lipidomics_data_modified[, 4:ncol(lipidomics_data_modified)], 
                  list(lipidomics_data_modified$Condition), 
                  mean) %>% 
  t(.) %>% 
  merge(x = ., y = lip_info, by = intersect(rownames(.), lip_info$Biochemical))
