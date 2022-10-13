# Download data from NCBI (Genome)
# Script by Rebecca E. Garner (rebecca.garner@mail.concordia.ca)

# Load libraries
library(tidyverse)
library(rentrez)


#### Import and format data ####
# Import (corrected) taxon list
tax_list <- read_tsv("output/2022-06-18_species_list.tsv", col_names = TRUE)


#### Download NCBI Genome data ####
# Get description of Genome database
entrez_db_summary("genome")

# List searchable fields for Genome database
entrez_db_searchable("genome")

# Download NCBI Genome data (https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html)
# Expected run time: 30-35 minutes
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
    ids <- entrez_search(db = "genome", term = paste0(tax_tmp[j], "[ORGN]"))$ids
    
    if (length(ids) > 0) {
      for (k in 1:length(ids)) {
        entrez_summary(db = "genome", id = ids[k]) %>%
          as_tibble() %>%
          mutate(species = tax_list$species[i],
                 query = tax_tmp[j]) %>%
          write_tsv(paste0("output/intermediate/ncbi_genome_records/ncbi_genome_records-", gsub(" ", "_", tax_list$species[i]), "-", j, "_id", ids[k], ".tsv"), col_names = TRUE, na = "")
      }
    }
  }
}
print(Sys.time() - start_time_tmp)

# Import and combine all NCBI genome records
filepath_ncbi_genome_records <- "output/intermediate/ncbi_genome_records/"

(ncbi_genome_records_files <- dir(path = filepath_ncbi_genome_records, pattern = "ncbi_genome_records-*"))
length(ncbi_genome_records_files)  # Number of files

ncbi_genome_records <- tibble(file_name = ncbi_genome_records_files) %>%
  mutate(file_contents = map(file_name, ~ read_tsv(file.path(filepath_ncbi_genome_records, .),
                                                   col_names = TRUE))) %>%
  unnest(c(file_contents))


#### Curate downloaded records ####
# Are any records missing query information?
any(is.na(ncbi_genome_records$species))  

# Any records empty?
any(is.na(ncbi_genome_records$uid))

# Are all NCBI Genome records for taxa in Desforges list?
tax_unique <- unique(c(tax_list$species, tax_list$synonym)) %>%
  str_split(", ") %>%
  unlist() %>%
  unique()

ncbi_genome_records %>%
  filter(!organism_name %in% tax_unique)

# Filter unique records
ncbi_genome_records <- ncbi_genome_records %>%
  distinct(species, uid, .keep_all = TRUE)

# Rearrange columns so species information is first
ncbi_genome_records <- ncbi_genome_records %>%
  relocate(c(species, query)) %>%
  relocate(file_name, .after = last_col())


#### Join conservation status information to BOLD records ####
# Join all columns in Desforges list OR BOLD records
bold_records_conservation <- tax_list %>%
  full_join(ncbi_genome_records, by = "species")

# Save NCBI genome record list to file.
# bold_records_conservation %>%
#   write_tsv(paste0("output/", Sys.Date(), "_ncbi_genome_records.tsv"), col_names = TRUE, na = "")
