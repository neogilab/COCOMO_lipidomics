##########################################################
# Script for cleaning data. This script will tidy the data 
# and handle missing values and zero or near-zero variance.
##########################################################

# Clear workspace
rm(list = ls())

# Load packages
library(tidyverse)
library(caret)

# Load data 
clinical_data <- read_csv("data/01_clinical_data.csv")
lipidomics_data <- read_csv("data/01_lipidomics_data.csv")
metabolomics_data <- read_csv("data/01_metabolomics_data.csv")
vatsat_data <- read_csv("data/01_vatsat_data.csv")



# Tidying clinical data  --------------------------------------------------
# Combine the "ID" col with the "Condition" col
clinical_data_2 <- as_tibble(clinical_data) %>%
  select(SUBJECT_ID, CONDITION, everything()) %>% 
  unite(col = "ID_Condition", "SUBJECT_ID":"CONDITION", sep = "_", remove = FALSE, na.rm = TRUE) 

# Join VAT/SAT data with the other clinical data
clinical_data_combined <- left_join(clinical_data_2, vatsat_data, by = "id")

# Subsetting the clinical variables
clinical_data_modified <- clinical_data_combined %>% 
  select(ID_Condition, AGE, GENDER, CONDITION, Ethnic, CD4_nadir, CDCAIDS, VAT, SAT,
         ART1, ART2, ART3, ART4, ART5, 
         ART1_prev, ART2_prev, ART3_prev, ART4_prev, ART5_prev, ART6_prev, ART7_prev, ART8_prev, ART9_prev, ART10_prev) %>% 
  replace_na(list(Ethnic = 4)) %>% # 4 = other/unknown
  replace_na(list(VAT = mean(.$VAT, na.rm = TRUE))) %>% 
  replace_na(list(SAT = mean(.$SAT, na.rm = TRUE))) %>% 
  #mutate(age_groupings = findInterval(AGE, c(40, 50, 60, 70))) %>% # Make age groupings 1[40-49], 2[50-59], 3[60-69], 4[70-79]
  mutate(Gender_enc = case_when(GENDER == 'Female' ~ 0,
                                GENDER == 'Male' ~ 1)) %>% 
  mutate(Condition_enc = case_when(CONDITION == 'HIV_NoMetS' ~ 0,
                                   CONDITION == 'HIV_MetS' ~ 1))

# Split the cART into their different groups
# Atripla
clinical_data_modified$ART1[clinical_data_modified$ART1 == "Atripla"] <- "Atripla_nrti"
clinical_data_modified$ART2[clinical_data_modified$ART1 == "Atripla_nrti"] <- "Atripla_nnrti"

# Genvoya
clinical_data_modified$ART1[clinical_data_modified$ART1 == "Genvoya"] <- "Genvoya_nrti"
clinical_data_modified$ART2[clinical_data_modified$ART1 == "Genvoya_nrti"] <- "Genvoya_pi"
clinical_data_modified$ART3[clinical_data_modified$ART1 == "Genvoya_pi"] <- "Genvoya_insti"

# Stribild
clinical_data_modified$ART1[clinical_data_modified$ART1 == "Stribild"] <- "Stribild_nrti"
clinical_data_modified$ART2[clinical_data_modified$ART1 == "Stribild_nrti"] <- "Stribild_pi"
clinical_data_modified$ART3[clinical_data_modified$ART1 == "Stribild_pi"] <- "Stribild_insti"

# Triumeq
clinical_data_modified$ART1[clinical_data_modified$ART1 == "Triumeq"] <- "Triumeq_nrti"
clinical_data_modified$ART2[clinical_data_modified$ART1 == "Triumeq_nrti"] <- "Triumeq_insti"

# List of unique names for all ART (curr + prev)
ART_all <- c(as_tibble(unique(clinical_data_modified$ART1))$value, as_tibble(unique(clinical_data_modified$ART2))$value, 
             as_tibble(unique(clinical_data_modified$ART3))$value, as_tibble(unique(clinical_data_modified$ART4))$value, 
             as_tibble(unique(clinical_data_modified$ART5))$value,
             as_tibble(unique(clinical_data_modified$ART1_prev))$value, as_tibble(unique(clinical_data_modified$ART2_prev))$value, 
             as_tibble(unique(clinical_data_modified$ART3_prev))$value, as_tibble(unique(clinical_data_modified$ART4_prev))$value, 
             as_tibble(unique(clinical_data_modified$ART5_prev))$value, as_tibble(unique(clinical_data_modified$ART6_prev))$value, 
             as_tibble(unique(clinical_data_modified$ART7_prev))$value, as_tibble(unique(clinical_data_modified$ART8_prev))$value, 
             as_tibble(unique(clinical_data_modified$ART9_prev))$value, as_tibble(unique(clinical_data_modified$ART10_prev))$value)
uniq_ART_all <- as_tibble(unique(ART_all))

# Based on the unique list of ART (43 terms), each term has been grouped into the following 6 groups
# The combined ART tablets: 'Atripla', 'Genvoya', 'Stribild', 'Triumeq' are split into the first 5 groups
NRTI            <- c('Combivir', 'Complera', 'Complera Complera', 'Descovy', 'Emtriva(emtricitabin,FTC)', 'Epivir(lamivudin,3TC)', 'Epzicom', 
                     'Kivexa', 'Lamivudin', 'Trizivir', 'Truvada', 'Viread(tenofovir,TDF)', 'Ziagen(abavacir,ABC)', 'Atripla_nrti', 'Genvoya_nrti', 'Stribild_nrti', 'Triumeq_nrti')
NNRTI           <- c('Edurant (rilpivirine, RPV)', 'Intelence (etravirine, ETR)', 'Invirase (saquinavir, SQV)', 'Stocrin', 
                     'Sustiva (Stocrin, efavirenz, EFV)', 'Viramune XR (nevirapine, NVP)', 'Atripla_nnrti')
PI              <- c('Amprenavir', 'Kaletra (Aluvia, lopinavir/ritonavir, LPV/r)', 'Norvir (ritonavir, RTV)', 'Norvir(ritonavir,RTV)','Prezcobix (Rezolsta, darunavir + cobicistat)', 
                     'Prezista (darunavir, DRV)', 'Reyataz (atazanavir, ATV)', 'Tybost(cobicstat)', 'Viracept (nelfinavir, NFV)', 'Genvoya_pi', 'Stribild_pi')
INSTI           <- c('Isentress (raltegravir)', 'Tivicay (dolutegravir)', 'Vitekta (elvitegravir)', 'Genvoya_insti', 'Stribild_insti', 'Triumeq_insti')
other_unknown   <- c('83', 'Fuzeon (enfuvirtide, ENF)', 'Selzentry (Celsentri, maraviroc)', 
                     'Retrovir(zidovudin,AZT)') # Only one lipo_art in "current ART cols", maybe a mistake - therefore in other group. Maraviroc is an entry inhibitor, however there are few of the, thus classified "other"

# List of previous ART
lipo_ART        <- c('Crixivan (indinavir, IDV)', 'Retrovir(zidovudin,AZT)', 'Videx(didanosine,ddl)', 'Videx(didanosine,ddl) Videx(didanosine)', 'Zerit(stavudin,d4T)')

# Creating new columns with encoded values for the ART groups and boolean column for thymidine exposure or not. 
druglist <- c('ART1', 'ART2', 'ART3', 'ART4', 'ART5')
clinical_data_modified_encoded <- clinical_data_modified %>% 
  mutate(NRTI = apply(.[, druglist], 1, function(r) any(r %in% NRTI)),
         NNRTI = apply(.[, druglist], 1, function(r) any(r %in% NNRTI)),
         PI = apply(.[, druglist], 1, function(r) any(r %in% PI)),
         INSTI = apply(.[, druglist], 1, function(r) any(r %in% INSTI)),
         Other_Unknown = apply(.[, druglist], 1, function(r) any(r %in% other_unknown)),
         thym_exposure = apply(.[, c('ART1_prev', 'ART2_prev', 'ART3_prev', 'ART4_prev', 'ART5_prev', 'ART6_prev', 'ART7_prev', 'ART8_prev', 'ART9_prev', 'ART10_prev')], 
                               1, function(r) any(r %in% lipo_ART)), # If a value from lipo_ART is in any of the 10 cols, then TRUE otherwise FALSE
         immunodeficiency = case_when(CD4_nadir > 200 | CDCAIDS == '0' ~ 0,
                                      CD4_nadir < 200 | CDCAIDS == '1' ~ 1)) %>% 
  select(ID_Condition, CONDITION, Condition_enc, AGE, GENDER, Gender_enc, Ethnic, immunodeficiency, thym_exposure, VAT, SAT,
         NRTI, NNRTI, PI, INSTI, Other_Unknown)

cols <- c('thym_exposure', 'NRTI', 'NNRTI', 'PI', 'INSTI', 'Other_Unknown')
clinical_data_modified_encoded[,cols] <- lapply(clinical_data_modified_encoded[,cols], as.numeric)



# Tidying lipidomics data ---------------------------------------------
# Remove rows as they contain clinical data, except subject_id
lipidomics_data_modified <- lipidomics_data %>% 
  slice(13:n()) 

# Save list with lipid info
lip_info <- lipidomics_data_modified %>% 
  select(BIOCHEMICAL, SUPER_PATHWAY,SUB_PATHWAY,	KEGG,	`Group   HMDB_ID`) %>% 
  slice(3:n()) %>% 
  rename(Biochemicals = BIOCHEMICAL, Super_pathway = SUPER_PATHWAY, Sub_pathway = SUB_PATHWAY, KEGG_ID = KEGG, HMDB_ID = `Group   HMDB_ID`) %>% 
  filter(!str_detect(Biochemicals, 'Total'))  %>% 
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
lipidomics_data_modified <- merge(x = lipidomics_data_modified, y = clinical_data_modified_encoded[ , c("ID_Condition", "GENDER")], by = "ID_Condition", all.x=TRUE) %>% 
#  filter(!GENDER == "Female") %>% 
  filter(!Condition == "Ctrl") %>% 
  select(GENDER, everything(), -ID)

# Make the chemical columns numeric
lipidomics_data_modified[, 4:ncol(lipidomics_data_modified)] <- sapply(lipidomics_data_modified[, 4:ncol(lipidomics_data_modified)], as.numeric)

# Check for missing values in lipid matrix
any(is.na(lipidomics_data_modified)) #FALSE = No missing values

# Remove columns with zero variance or near zero variance
compute_nzv <- nearZeroVar(lipidomics_data_modified[, 4:ncol(lipidomics_data_modified)], 
                           saveMetrics = TRUE) %>% 
  mutate(Biochemicals = rownames(.))

nzv_lipids <- compute_nzv %>% 
  filter(nzv == FALSE) %>% 
  select(Biochemicals)

lipidomics_data_modified_var <- lipidomics_data_modified %>% 
  select(GENDER, ID_Condition, Condition, all_of(nzv_lipids$Biochemicals))

# Change TAG names to accomodate the correct nomenclature (insert paranthesis)
col_names <- gsub('^(TAG)([0-9]*\\S[0-9]*)(\\SFA)([0-9]*\\S[0-9]*)$', '\\1(\\2)\\3(\\4)', colnames(lipidomics_data_modified_var))
colnames(lipidomics_data_modified_var) <- col_names

# Change TAG names
lip_info_without_nzv <- lip_info %>% 
  filter(Biochemicals %in% nzv_lipids$Biochemicals) %>% 
  mutate(Biochemicals = col_names[4:920])



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
  rename(Biochemicals = ...2, KEGG_ID = ...12, HMDB_ID = `CLIENT IDENTIFIER`, Super_pathway = ...3, Sub_pathway = ...4) %>% 
  filter(Biochemicals %in% sign_met)  %>% 
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

# Make the metabolite columns numeric
metabolomics_data_modified[, 4:ncol(metabolomics_data_modified)] <- sapply(metabolomics_data_modified[, 4:ncol(metabolomics_data_modified)], as.numeric)

# Check for missing values in the metabolite matrix
any(is.na(metabolomics_data_modified)) #FALSE = No missing values

# Remove columns with zero variance or near zero variance
compute_nzv <- nearZeroVar(metabolomics_data_modified[, 4:ncol(metabolomics_data_modified)], saveMetrics = TRUE) %>% 
  mutate(Biochemicals = rownames(.))

nzv_metabolites <- compute_nzv %>% 
  filter(nzv == FALSE) %>% 
  select(Biochemicals)

metabolomics_data_modified_var <- metabolomics_data_modified %>% 
  select(ID_Condition, Condition, nzv_metabolites$Biochemicals)

# Select significant metabolites
metabolomics_sign_met <- metabolomics_data_modified %>% 
  select(ID_Condition, sign_met)

  

# Save data to csv  -------------------------------------------
lipidomics_data_modified_var %>% 
  write_csv("data/02_lipidomics_data_tidy.csv")

lip_info_without_nzv %>% 
  write_csv("data/02_lipidomics_data_info.csv")

metabolomics_data_modified_var %>% 
  write_csv("data/02_metabolomics_data_tidy.csv")

met_info %>% 
  write_csv("data/02_metabolomics_data_info.csv")

metabolomics_sign_met %>% 
  write_csv("data/02_metabolomics_sign_metabolites.csv")

clinical_data_modified_encoded %>% 
  write_csv("data/02_clinical_data_tidy.csv")
