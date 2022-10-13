# Download data from BOLD Systems
# Script by Rebecca E. Garner (rebecca.garner@mail.concordia.ca)

# Load libraries
library(tidyverse)
library(bold)


#### Import and format data ####
# Import (corrected) taxon list
tax_list <- read_tsv("output/2022-06-18_species_list.tsv", col_names = TRUE)


#### Download BOLD records ####
# Download specimen and sequence data (https://docs.ropensci.org/bold/)
start_time_tmp <- Sys.time()
for (i in 1:nrow(tax_list)) {
  print(paste0("Now processing taxon number ", i, " (", tax_list$species[i], ")"))
  
  if (is.na(tax_list$synonym[i])) {
    tax_tmp <- tax_list$species[i]
  } else if (!is.na(tax_list$synonym[i])) {
    tax_tmp <- unique(c(tax_list$species[i], tax_list$synonym[i])) %>%
      str_split(", ") %>%
      unlist() %>%
      unique()
  }
  
  for (j in 1:length(tax_tmp)) {
    bold_seqspec(taxon = tax_tmp[j]) %>%
      as_tibble() %>%
      mutate(species = tax_list$species[i],
             query = tax_tmp[j]) %>%
      write_tsv(paste0("output/intermediate/bold_records/bold_records-", gsub(" ", "_", tax_list$species[i]), "-", j, ".tsv"), col_names = TRUE, na = "")
  }
}
print(Sys.time() - start_time_tmp)

# Import and combine all bold records files
filepath_bold_records <- "output/intermediate/bold_records/"

(bold_records_files <- dir(path = filepath_bold_records, pattern = "bold_records-*"))
length(bold_records_files)  # Number of files

bold_records <- tibble(file_name = bold_records_files) %>%
  mutate(file_contents = map(file_name, ~ read_tsv(file.path(filepath_bold_records, .),
                                                   col_names = TRUE, col_types = cols(.default = "c")))) %>%
  unnest(c(file_contents))

# Were records not downloaded for any taxa?
unique(tax_list$species)[which(!unique(tax_list$species) %in% unique(bold_records$species))]

# Import manually downloaded records
filepath_bold_records_manual <- "output/intermediate/bold_records_manual/"
(bold_records_manual_files <- dir(path = filepath_bold_records_manual, pattern = "^bold_data-*"))
length(bold_records_manual_files)  # Number of files

bold_records_manual <- tibble(file_name = bold_records_manual_files) %>%
  mutate(file_contents = map(file_name, ~ read_tsv(file.path(filepath_bold_records_manual, .),
                                                   col_names = TRUE, col_types = cols(.default = "c")))) %>%
  unnest(c(file_contents))

bold_records_manual$species <- NA
bold_records_manual$query <- NA
for (i in 1:nrow(bold_records_manual)) {
  tax_tmp <- bold_records_manual$file_name[i] %>%
    str_remove("^bold_data-") %>%
    str_remove("-.*") %>%
    str_replace("_", " ")
  
  query_n_tmp <- bold_records_manual$file_name[i] %>%
    str_remove("^bold_data-") %>%
    str_remove(".txt$") %>%
    str_remove(".*-") %>%
    as.numeric()
  
  if (query_n_tmp == 1) {
    bold_records_manual$species[i] <- tax_list$species[which(tax_list$species == tax_tmp)]
    bold_records_manual$query[i] <- tax_list$species[which(tax_list$species == tax_tmp)]
  } else if (query_n_tmp > 1) {
    bold_records_manual$species[i] <- tax_list$species[which(tax_list$species == tax_tmp)]
    
    synonyms_tmp <- tax_list$synonym[which(tax_list$species == tax_tmp)] %>%
      str_split(", ") %>%
      unlist()
    
    bold_records_manual$query[i] <- synonyms_tmp[query_n_tmp - 1]
  }
}

# Collate programmatically and manually downloaded records
bold_records_all <- bold_records %>%
  filter(!query %in% bold_records_manual$query) %>%
  bind_rows(bold_records_manual)


#### Curate downloaded records ####
# Are any records missing query information?
any(is.na(bold_records_all$species))  

if (any(is.na(bold_records_all$species))) {
  for (i in 1:nrow(bold_records_all)) {
    if (is.na(bold_records_all$species[i]) | is.na(bold_records_all$query[i])) {
      filename_tmp <- gsub("^bold_data-", "bold_records-", bold_records_all$file_name[i])
      filename_tmp <- gsub(".txt$", ".tsv", filename_tmp)
      
      tax_tmp <- filename_tmp %>%
        str_remove("^bold_records-") %>%
        str_remove("-.*") %>%
        str_replace("_", " ")
      
      query_n_tmp <- filename_tmp %>%
        str_remove("^bold_records-") %>%
        str_remove(".tsv$") %>%
        str_remove(".*-") %>%
        as.numeric()
      
      if (query_n_tmp == 1) {
        bold_records_all$species[i] <- tax_list$species[which(tax_list$species == tax_tmp)]
        bold_records_all$query[i] <- tax_list$species[which(tax_list$species == tax_tmp)]
      } else if (query_n_tmp > 1) {
        bold_records_all$species[i] <- tax_list$species[which(tax_list$species == tax_tmp)]
        
        synonyms_tmp <- tax_list$synonym[which(tax_list$species == tax_tmp)] %>%
          str_split(", ") %>%
          unlist()
        
        bold_records_all$query[i] <- synonyms_tmp[query_n_tmp - 1]
      }
    }
  }
}

# Delete empty records (processid is unique identifier in BOLD)
bold_records_all <- bold_records_all %>%
  filter(!is.na(processid))

# Filter unique records
bold_records_all <- bold_records_all %>%
  distinct(species, processid, .keep_all = TRUE)
length(unique(bold_records_all$processid))  # Number of unique records

# Inspect BOLD process IDs duplicated in records list
nrow(bold_records_all) - length(unique(bold_records_all$processid))

bold_processid_dup <- bold_records_all %>%
  group_by(processid) %>%
  dplyr::count() %>%
  filter(n > 1) %>%
  pull(processid)

bold_records_all %>%
  filter(processid %in% bold_processid_dup) %>%
  dplyr::select(processid, species, query) %>%
  arrange(processid)

# Rearrange columns so species information is first
bold_records_all <- bold_records_all %>%
  relocate(c(species, query)) %>%
  relocate(file_name, .after = last_col()) %>%
  dplyr::select(-value)


#### Join conservation status information to BOLD records ####
# Join all columns in Desforges list OR BOLD records
bold_records_conservation <- tax_list %>%
  full_join(bold_records_all, by = "species")

# Are any taxa not represented in records list?
unique(tax_list$species)[which(!unique(tax_list$species) %in% unique(bold_records_conservation$species))]

# Save BOLD records with conservation status information to file
# bold_records_conservation %>%
#   write_tsv(paste0("output/", Sys.Date(), "_bold_records.tsv"), col_names = TRUE, na = "")
