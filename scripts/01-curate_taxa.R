# Curate taxon list
# Script by Rebecca E. Garner (rebecca.garner@mail.concordia.ca)

# Load libraries
library(tidyverse)
library(janitor)
library(taxize)
library(natserv)


#### Import and format data ####
# Download Desforges et al. (2022) species list from
# https://cdnsciencepub.com/doi/suppl/10.1139/cjfas-2021-0073/suppl_file/cjfas-2021-0073suppla.csv

# Import and format Desforges list, removing weird characters, whitespace, and
# the word "applicable" which shows up in a few taxon names
spp_list <- read_csv("data/cjfas-2021-0073suppla.csv", skip = 46, col_names = TRUE) %>%
  clean_names() %>%
  rename(taxon = taxa,
         assessment = overall_assessment) %>%
  mutate(across(everything(), ~str_replace(.x, "\u00A0", " "))) %>%  # Replace no-breaking space with space
  mutate(across(everything(), ~str_remove_all(.x, "\u00C2"))) %>%  # Remove capital A with circumflex accent
  mutate(across(everything(), ~str_remove_all(.x, "\u00AC"))) %>%  # Remove not sign
  mutate(across(everything(), ~str_remove_all(.x, "\u2020"))) %>%  # Remove dagger
  mutate(across(everything(), ~str_trim(.x, "both"))) %>%  # Trim leading and trailing whitespace
  mutate(across(c(taxon, order, family, species), ~str_to_sentence(.x))) %>%
  mutate(source = str_to_upper(source)) %>%
  mutate(assessment = str_to_title(assessment)) %>%
  mutate(species = str_remove(species, " applicable"))


#### Resolve duplicate species entries ####
# Identify species with duplicate entries in Desforges list
spp_duplicates <- spp_list %>%
  group_by(species) %>%
  count() %>%
  filter(n > 1) %>%
  pull(species)

# Resolve duplicate or conflicting conservation status assessments
spp_list_duplicates_collapsed <- spp_list %>%
  filter(species %in% spp_duplicates) %>%
  distinct() %>%
  group_by(taxon, order, family, species, origin, assessed, source) %>%
  summarize(assessment = toString(assessment)) %>%
  ungroup() %>%
  mutate(assessment = case_when(!grepl(", ", assessment) ~ assessment,
                                grepl("Endangered", assessment) & grepl("Special Concern", assessment) ~ "Endangered",
                                grepl("Endangered", assessment) & grepl("Threatened", assessment) ~ "Endangered",
                                TRUE ~ assessment))

# Update species lists with resolved duplicated species
spp_list <- spp_list %>%
  filter(!species %in% spp_duplicates) %>%
  bind_rows(spp_list_duplicates_collapsed)


#### Look up taxon spellings in databases #####
# List data sources accessible to taxize (under the column "title")
gnr_datasources()

# Extract database ID(s)
refdb <- subset(gnr_datasources(), title %in% c("Integrated Taxonomic Information SystemITIS",
                                                "GBIF Backbone Taxonomy")) %>%
  pull(id)

# Compare spelling in Desforges list and databases
# Expected run time: <10 minutes
start_time_tmp <- Sys.time()
tax_spelling_all <- tibble(NULL)
for (i in 1:nrow(spp_list)) {
  tax_spelling_tmp <- spp_list$species[i] %>%
    gnr_resolve(data_source_ids = refdb, with_canonical_ranks = TRUE)
  tax_spelling_all <- tax_spelling_all %>%
    bind_rows(tax_spelling_tmp)
}
print(Sys.time() - start_time_tmp)

# Write taxon name matching table to file
# tax_spelling_all %>%
#   write_tsv(paste0("output/intermediate/", Sys.Date(), "_taxize_spelling_", length(refdb), "dbs.tsv"), col_names = TRUE)


#### Assess taxon spellings ####
# Import taxon name matching table from file
tax_spelling_all <- read_tsv("output/intermediate/2022-06-16_taxize_spelling_2dbs.tsv", col_names = TRUE)

# Assign taxonomic ranks to matches
tax_spelling_all <- tax_spelling_all %>%
  mutate(match_tax_rank = case_when(str_count(matched_name2, " ") == 0 ~ "genus",
                                    str_count(matched_name2, " ") == 1 ~ "species",
                                    str_count(matched_name2, " ") == 2 ~ "subspecies"),
         desforges_tax_rank = case_when(str_count(user_supplied_name, " ") == 0 ~ "genus",
                                        str_count(user_supplied_name, " ") == 1 ~ "species",
                                        str_count(user_supplied_name, " ") == 2 ~ "subspecies"))

# If match only returns a genus name, does the genus match the genus in Desforges?
tax_spelling_all %>%
  filter(match_tax_rank == "genus") %>%
  mutate(matching_genus = case_when(str_remove_all(user_supplied_name, " .*") == matched_name2 ~ "yes",
                                    TRUE ~ "no")) %>%
  filter(matching_genus == "no")  # If returns empty tibble, then all genera match

# If only returns a species name when a subspecies is supplied, does the species match the species in Desforges?
tax_spelling_all %>%
  filter(match_tax_rank == "species" & desforges_tax_rank == "subspecies") %>%
  separate(user_supplied_name, into = c("genus_desforges", "species_desforges", "subspecies_desforges"), sep = " ") %>%
  mutate(matching_species = case_when(str_c(genus_desforges, species_desforges, sep = " ") == matched_name2 ~ "yes",
                                      TRUE ~ "no")) %>%
  filter(matching_species == "no")  # If returns empty tibble, then all species match

# Summarize taxon spelling matches
tax_spelling_matches <- tax_spelling_all %>%
  distinct(user_supplied_name, score, matched_name2, match_tax_rank,desforges_tax_rank) %>%
  mutate(match_agreement = case_when(user_supplied_name == matched_name2 ~ "match",
                                     user_supplied_name != matched_name2 & match_tax_rank == "genus" & str_remove_all(user_supplied_name, " .*") == matched_name2 ~ "mismatch_genus",
                                     user_supplied_name != matched_name2 & match_tax_rank == "species" & desforges_tax_rank == "species" ~ "mismatch_species",
                                     user_supplied_name != matched_name2 & desforges_tax_rank == "subspecies" ~ "mismatch_subspecies",
                                     TRUE ~ "mismatch_other"))

# How many taxa have multiple returned matches from taxize?
tax_spelling_all %>%
  distinct(user_supplied_name, score, matched_name2) %>%
  group_by(user_supplied_name) %>%
  count() %>%
  filter(n > 1)

# Do some entries have different types of matches from different databases?
tax_spelling_matches %>%
  distinct(user_supplied_name, match_agreement) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = match_agreement, values_from = value, values_fill = 0) %>%
  mutate(match_types = match + mismatch_species + mismatch_genus + mismatch_subspecies) %>%
  filter(match_types > 1)

# Remove matches for genera only when an exact match is also returned
tax_match_exact_plus_genus <- tax_spelling_matches %>%
  distinct(user_supplied_name, match_agreement) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = match_agreement, values_from = value, values_fill = 0) %>%
  mutate(match_types = match + mismatch_species + mismatch_genus + mismatch_subspecies) %>%
  filter(match >= 1 & mismatch_genus >= 1) %>%
  pull(user_supplied_name)

tax_spelling_matches <- tax_spelling_matches %>%
  filter(!(match_tax_rank == "genus" & user_supplied_name %in% tax_match_exact_plus_genus))

# Extract taxa not returned by Taxize
no_matches <- spp_list %>%
  filter(!species %in% unique(tax_spelling_all$user_supplied_name)) %>%
  pull(species)
length(no_matches)  # Number of taxa not returned by Taxize

# Extract taxa with matched spellings
matched_spellings <- tax_spelling_matches %>%
  filter(match_agreement == "match") %>%
  distinct(user_supplied_name, matched_name2) %>%
  pull(user_supplied_name)
length(matched_spellings)  # Number of taxa with matched spellings

# Extract taxa with mismatched genus
mismatch_genus_spellings <- tax_spelling_matches %>%
  filter(match_agreement == "mismatch_genus") %>%
  distinct(user_supplied_name, matched_name2) %>%
  pull(user_supplied_name)
length(mismatch_genus_spellings)  # Number of taxa with mismatched genus

# Extract taxa with mismatched subspecies
mismatch_subspecies_spellings <- tax_spelling_matches %>%
  filter(match_agreement == "mismatch_subspecies") %>%
  distinct(user_supplied_name, matched_name2) %>%
  pull(user_supplied_name)
length(mismatch_subspecies_spellings)  # Number of taxa with mismatched subspecies

# Summarize alternate spellings for mismatches
alt_spellings <- tax_spelling_matches %>%
  filter(match_agreement == "mismatch_species") %>%
  group_by(user_supplied_name) %>%
  summarize(taxize_spelling = toString(matched_name2)) %>%
  ungroup()
length(alt_spellings$user_supplied_name)  # Number of taxa with alternate spellings from matches

# Join Desforges list with Taxize information
tax_spelling_summary <- spp_list %>%
  mutate(taxize = case_when(species %in% no_matches ~ "no_match",
                            species %in% matched_spellings ~ "matched",
                            species %in% mismatch_genus_spellings ~ "mismatch_genus",
                            species %in% mismatch_subspecies_spellings ~ "mismatch_subspecies",
                            species %in% alt_spellings$user_supplied_name ~ "alternate_spelling")) %>%
  left_join(alt_spellings, by = c("species" = "user_supplied_name")) %>%
  rename(species_desforges = species) %>%
  arrange(taxon, order, family, species_desforges)

# Write table of unresolved taxon spellings
# tax_spelling_summary %>%
#   filter(taxize != "matched") %>%
#   write_csv(paste0("output/intermediate/", Sys.Date(), "_unresolved_spellings.csv"), col_names = TRUE, na = "")


#### Manually curate misspellings ####
# Manually correct taxon misspellings
tax_curated <- tax_spelling_summary %>%
  mutate(species_corrected = case_when(grepl("^Acerpenot ", species_desforges) ~ str_replace(species_desforges, "^Acerpenot ", "Acerpenna "),
                                       species_desforges == "Amystoma maculatum" ~ "Ambystoma maculatum",
                                       species_desforges == "Asynarchus cinot moneus" ~ "Asynarchus cinnamoneus",
                                       species_desforges == "Enallagma anot" ~ "Enallagma anna",
                                       species_desforges == "Enallagma antenot tum" ~ "Enallagma antennatum",
                                       species_desforges == "Esatina eschscholtzii" ~ "Ensatina eschscholtzii",
                                       species_desforges == "Fisheserola nuttallii" ~ "Fisherola nuttalli",
                                       species_desforges == "Graptemy geographica" ~ "Graptemys geographica",
                                       species_desforges == "Hyla chysoscelis" ~ "Hyla chrysoscelis",
                                       species_desforges == "Lithobates sylvacticus" ~ "Lithobates sylvaticus",
                                       grepl("^Molanot ", species_desforges) ~ str_replace(species_desforges, "^Molanot ", "Molanna "),
                                       species_desforges == "Plegais falcinellus" ~ "Plegadis falcinellus",
                                       species_desforges == "Wormaldia clauseni" ~ "Wormaldia occidea clauseni",
                                       TRUE ~ species_desforges)) %>%
  filter(!species_desforges %in% c("Olivaria subrotunda"))


#### Expand list of taxonomic synonyms ####
# Add column for taxonomic synonyms listed in NatureServe Explorer
# Expected run time: ~20 minutes
start_time_tmp <- Sys.time()
tax_curated$synonym <- NA
for (i in 1:nrow(tax_curated)) {
  species_corrected_tmp <- str_split(tax_curated$species_corrected[i], ", ")[[1]]
  print(paste0("Now processing taxon number ", i, " (", toString(species_corrected_tmp), ")"))
  
  species_all_tmp <- sub("^(\\S*\\s+\\S+).*", "\\1", species_corrected_tmp)
  
  for (j in 1:length(species_corrected_tmp)) {
    synonym_tmp <- ns_search_spp(species_corrected_tmp[j])$results$speciesGlobal$synonyms[[1]]
    synonym_tmp <- gsub("\\sx\\s", " ", synonym_tmp)
    synonym_tmp <- sub("^(\\S*\\s+\\S+).*", "\\1", synonym_tmp)
    
    if (length(synonym_tmp) >= 1) {
      for (k in 1:length(synonym_tmp)) {
        if (!synonym_tmp[k] %in% species_all_tmp &
            !synonym_tmp[k] %in% unique(tax_curated$species_corrected) &
            synonym_tmp[k] != "" &
            !grepl("\\w+\\ssp\\.", synonym_tmp[k])) {
          species_all_tmp <- c(species_all_tmp, synonym_tmp[k])
        }
      }
    }
  }
  tax_curated$synonym[i] <- species_all_tmp[species_all_tmp != species_corrected_tmp & !species_all_tmp %in% unique(tax_curated$species_corrected)] %>%
    toString()
}
tax_curated$synonym[tax_curated$synonym == ""] <- NA
print(Sys.time() - start_time_tmp)

# Write taxonomic synonym table to file
# tax_curated %>%
#   write_tsv(paste0("output/intermediate/", Sys.Date(), "_tax_synonyms.tsv"), col_names = TRUE, na = "")


#### Rewrite taxonomy upstream of genus ####
tax_curated <- read_tsv("output/intermediate/2022-06-16_tax_synonyms.tsv", col_names = TRUE)

# Create vector of unique genus+species
desforges_species <- unique(sub("^(\\S*\\s+\\S+).*", "\\1", unique(tax_curated$species_corrected)))

# Assign taxonomy with Taxize classification
# **WARNING: Requires user input when taxa return multiple TSNs**
start_time_tmp <- Sys.time()
tax_ranks <- tibble(NULL)
for (i in 1:length(desforges_species)) {
  print(paste0("Now processing taxon number ", i, " (", desforges_species[i], ")"))
  
  species_tmp <- desforges_species[i]
  
  classification_tmp <- classification(species_tmp, db = "itis", accepted = TRUE, ask = TRUE)[[1]]
  
  if (any(!is.na(classification_tmp))) {
    classification_tmp <- classification_tmp %>%
      dplyr::select(name, rank) %>%
      pivot_wider(names_from = rank, values_from = name)
    
    genus_classification_tmp <- bind_cols(query = species_tmp, classification_tmp) %>%
      mutate(db = "ITIS")
    
    tax_ranks <- tax_ranks %>%
      bind_rows(genus_classification_tmp)
  } else {
    classification_tmp <- classification(species_tmp, db = "gbif", accepted = TRUE, ask = TRUE)[[1]]
    
    if (any(!is.na(classification_tmp))) {
      classification_tmp <- classification_tmp %>%
        dplyr::select(name, rank) %>%
        pivot_wider(names_from = rank, values_from = name)
      
      genus_classification_tmp <- bind_cols(query = species_tmp, classification_tmp) %>%
        mutate(db = "GBIF")
      
      tax_ranks <- tax_ranks %>%
        bind_rows(genus_classification_tmp)
    }
  }
}
print(Sys.time() - start_time_tmp)

# Write taxonomy to file
# tax_ranks %>%
#   write_tsv(paste0("output/intermediate/", Sys.Date(), "_classifications.tsv"), col_names = TRUE, na = "")

tax_ranks <- read_tsv("output/intermediate/2022-06-16_classifications.tsv", col_names = TRUE, col_types = cols(.default = "c"))

# Summarize taxonomy by database
tax_ranks_db <- tax_ranks %>%
  filter(!(query == "Chara aspera" & species == "Bolitochara aspera")) %>%  # Remove apparently erroneous taxonomy
  mutate(genus_query = str_remove_all(query, " .*")) %>%
  dplyr::select(query, genus_query, db, class, order, family, genus, species) %>%
  distinct(genus_query, db, class, order, family) %>%
  pivot_longer(!c(genus_query, db), names_to = "rank", values_to = "taxonomy") %>%
  mutate(db = str_to_lower(db)) %>%
  pivot_wider(names_from = c(rank, db), values_from = taxonomy)

# Verify completeness of database taxonomy
tax_ranks_db$genus_query[which(is.na(tax_ranks_db$class_itis) & is.na(tax_ranks_db$class_gbif))]
tax_ranks_db$genus_query[which(is.na(tax_ranks_db$order_itis) & is.na(tax_ranks_db$order_gbif))]
tax_ranks_db$genus_query[which(is.na(tax_ranks_db$family_itis) & is.na(tax_ranks_db$family_gbif))]

# Manually fill missing taxa
tax_ranks_db <- tax_ranks_db %>%
  mutate(order_itis = case_when(genus_query == "Haitia" ~ "Basommatophora",
                                genus_query == "Ladislavella" ~ "Basommatophora",
                                genus_query == "Sibirenauta" ~ "Basommatophora",
                                TRUE ~ order_itis))

# Replace Desforges taxonomy upstream of genus
tax_curated$genus <- NA
tax_curated$family_curated <- NA
tax_curated$order_curated <- NA
tax_curated$class_curated <- NA
for (i in 1:nrow(tax_curated)) {
  tax_curated$genus[i] <- str_remove_all(tax_curated$species_corrected[i], " .*")
  
  family_tmp <- tax_ranks_db$family_itis[which(tax_ranks_db$genus_query == tax_curated$genus[i])]
  if (is.na(family_tmp)) {
    family_tmp <- tax_ranks_db$family_gbif[which(tax_ranks_db$genus_query == tax_curated$genus[i])]
  }
  tax_curated$family_curated[i] <- family_tmp
  
  order_tmp <- tax_ranks_db$order_itis[which(tax_ranks_db$genus_query == tax_curated$genus[i])]
  if (is.na(order_tmp)) {
    order_tmp <- tax_ranks_db$order_gbif[which(tax_ranks_db$genus_query == tax_curated$genus[i])]
  }
  tax_curated$order_curated[i] <- order_tmp
  
  class_tmp <- tax_ranks_db$class_itis[which(tax_ranks_db$genus_query == tax_curated$genus[i])]
  if (is.na(class_tmp)) {
    class_tmp <- tax_ranks_db$class_gbif[which(tax_ranks_db$genus_query == tax_curated$genus[i])]
  }
  tax_curated$class_curated[i] <- class_tmp
}

# Inspect upstream taxonomic agreement
tax_curated %>%
  distinct(genus, family_curated) %>%
  group_by(genus) %>%
  count() %>%
  filter(n > 1)

tax_curated %>%
  distinct(family_curated, order_curated) %>%
  group_by(family_curated) %>%
  count() %>%
  filter(n > 1)  # Families with disagreement in order-rank taxonomy

tax_curated %>%
  distinct(order_curated, class_curated) %>%
  group_by(order_curated) %>%
  count() %>%
  filter(n > 1)  # Orders with disagreement in class-rank taxonomy

tax_curated %>%
  distinct(class_curated, taxon) %>%
  group_by(class_curated) %>%
  count() %>%
  filter(n > 1)

# Curate and harmonize database-informed taxonomy
tax_curated <- tax_curated %>%
  mutate(order_curated = case_when(family_curated == "Amnicolidae" ~ "Neotaenioglossa",
                                   family_curated == "Unionidae" ~ "Unionoida",
                                   TRUE ~ order_curated),
         class_curated = case_when(order_curated == "Alismatales" ~ "Magnoliopsida",
                                   order_curated == "Asparagales" ~ "Magnoliopsida",
                                   order_curated == "Poales" ~ "Magnoliopsida",
                                   TRUE ~ class_curated))

# Write curated taxonomy to file
# tax_curated %>%
#   write_tsv(paste0("output/intermediate/", Sys.Date(), "_taxonomy.tsv"), col_names = TRUE, na = "")


#### Write finalized curated species list to file ####
tax_curated <- read_tsv("output/intermediate/2022-06-18_taxonomy.tsv", col_names = TRUE)

# Write spelling-curated taxon list to file
# tax_curated %>%
#   select(taxon, class_curated, order_curated, family_curated, genus,
#          species_corrected, synonym, origin, assessed, source, assessment) %>%
#   rename(class = class_curated,
#          order = order_curated,
#          family = family_curated,
#          species = species_corrected) %>%
#   arrange(taxon, class, order, family, genus, species) %>%
#   write_tsv(paste0("output/", Sys.Date(), "_species_list.tsv"), col_names = TRUE, na = "")
