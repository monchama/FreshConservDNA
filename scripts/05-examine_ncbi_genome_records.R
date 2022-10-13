# Examine NCBI genome records
# Script by Rebecca E. Garner (rebecca.garner@mail.concordia.ca)

# Load libraries
library(tidyverse)

# Load palettes
source("scripts/00_palettes.R")

#### Import and format data ####
# Import BOLD records
ncbigenome_records <- read_tsv("output/2022-06-18_ncbi_genome_records.tsv", col_names = TRUE, col_types = cols(.default = "c"))

# Extract taxonomy and conservation status information
conservation_status <- ncbigenome_records %>%
  distinct(taxon, class, order, family, genus, species, assessment)


#### Count taxa in Desforges list ####
# Determine number of taxa by "taxon" (Amphibians, Birds, Fishes, ...) in Desforges list
(n_taxa_by_taxon_desforges <- conservation_status %>%
   group_by(taxon) %>%
   count(name = "n_taxa_desforges") %>%
   ungroup() %>%
   arrange(-n_taxa_desforges))

# Determine number of taxa by order in Desforges list
(n_taxa_by_order_desforges <- conservation_status %>%
    group_by(taxon, class, order) %>%
    count(name = "n_taxa_desforges") %>%
    ungroup() %>%
    arrange(-n_taxa_desforges))

# Determine number of taxa by conservation status
(n_taxa_by_assessment_desforges <- conservation_status %>%
    group_by(assessment) %>%
    count(name = "n_taxa_desforges") %>%
    ungroup() %>%
    arrange(-n_taxa_desforges))


#### Examine conservation statuses of freshwater taxa ####
# Number of species per conservation status category
(n_species_conservation <- conservation_status %>%
   group_by(taxon, assessment) %>%
   count(name = "n_species") %>%
   ungroup() %>%
   pivot_wider(names_from = taxon, values_from = n_species, values_fill = 0) %>%
   arrange(factor(assessment, levels = cosewic_ranking)) %>%
   mutate(n_species_total = rowSums(across(where(is.numeric)))))


#### Examine NCBI genome records ####
# Determine number of total vs. unique NCBI genome records
(n_ncbigenome_records_total <- ncbigenome_records %>%
   filter(!is.na(uid)) %>%
   nrow())  # Number of total NCBI genome records

(n_ncbigenome_records_unique <- ncbigenome_records %>%
    filter(!is.na(uid)) %>%
    distinct(uid, .keep_all = TRUE) %>%
    nrow())  # Number of unique NCBI genome records

# Determine species without records
(tax_without_records <- ncbigenome_records %>%
    filter(is.na(uid)) %>%
    distinct(species) %>%
    pull(species))
length(tax_without_records)  # Number of species without records

# Determine species with records
taxa_with_records <- ncbigenome_records %>%
  filter(!is.na(uid)) %>%
  distinct(species) %>%
  pull(species)
length(taxa_with_records)  # Number of species with records

# Count total NCBI genome records by species
n_ncbigenome_records <- ncbigenome_records %>%
  mutate(record = case_when(is.na(uid) ~ 0,
                            TRUE ~ 1)) %>%
  group_by(taxon, class, order, family, genus, species) %>%
  summarize(n_ncbigenome_records = sum(record),
            ncbigenome_assembly_accession = toString(assembly_accession)) %>%
  ungroup()


#### Examine genome assembly status ####
# Determine genome assembly statuses across NCBI genome records
ncbigenome_records %>%
  filter(!is.na(uid)) %>%
  distinct(uid, .keep_all = TRUE) %>%
  group_by(status) %>%
  count(name = "n_records")


#### Examine NCBI genome records by taxonomic group ####
# Determine number of taxa with NCBI genome records by "taxon"
(n_taxa_by_taxon_ncbigenome <- ncbigenome_records %>%
   filter(!is.na(uid)) %>%
   distinct(uid, .keep_all = TRUE) %>%
   distinct(taxon, species) %>%
   group_by(taxon) %>%
   count(name = "n_taxa_ncbigenome") %>%
   ungroup() %>%
   arrange(-n_taxa_ncbigenome))
sum(n_taxa_by_taxon_ncbigenome$n_taxa_ncbigenome) == n_ncbigenome_records_unique  # Should evaluate to TRUE

# Determine number of taxa with NCBI genome records by order
(n_taxa_by_order_ncbigenome <- ncbigenome_records %>%
    filter(!is.na(uid)) %>%
    distinct(uid, .keep_all = TRUE) %>%
    distinct(taxon, class, order, species) %>%
    group_by(taxon, class, order) %>%
    count(name = "n_taxa_ncbigenome") %>%
    ungroup() %>%
    arrange(-n_taxa_ncbigenome))
sum(n_taxa_by_order_ncbigenome$n_taxa_ncbigenome) == n_ncbigenome_records_unique  # Should evaluate to TRUE

# Summarize NCBI genome records by "taxon"
(ncbigenome_records_summary_by_taxon <- ncbigenome_records %>%
    filter(!is.na(uid)) %>%
    distinct(uid, .keep_all = TRUE) %>%
    group_by(taxon) %>%
    count(name = "n_records") %>%
    ungroup() %>%
    arrange(-n_records) %>%
    left_join(n_taxa_by_taxon_ncbigenome, "taxon") %>%
    left_join(n_taxa_by_taxon_desforges, "taxon") %>%
    mutate(pct_taxa_ncbigenome = n_taxa_ncbigenome/n_taxa_desforges * 100))

# Summarize NCBI genome records by order
(ncbigenome_records_summary_by_order <- ncbigenome_records %>%
    filter(!is.na(uid)) %>%
    distinct(uid, .keep_all = TRUE) %>%
    group_by(taxon, class, order) %>%
    count(name = "n_records") %>%
    ungroup() %>%
    arrange(-n_records) %>%
    left_join(n_taxa_by_order_ncbigenome, c("taxon", "class", "order")) %>%
    left_join(n_taxa_by_order_desforges, c("taxon", "class", "order")) %>%
    mutate(pct_taxa_ncbigenome = n_taxa_ncbigenome/n_taxa_desforges * 100))


#### Examine NCBI genome records by conservation status ####
# Determine NCBI genome records by conservation status assessment
(ncbigenome_records_summary_by_assessment <- ncbigenome_records %>%
   filter(!is.na(uid)) %>%
   group_by(assessment) %>%
   count(name = "n_records") %>%
   ungroup() %>%
   left_join(n_taxa_by_assessment_desforges, "assessment") %>%
   mutate(pct_records_ncbigenome = n_records/n_taxa_desforges * 100) %>%
   arrange(-pct_records_ncbigenome))
