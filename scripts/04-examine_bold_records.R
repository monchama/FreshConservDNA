# Examine BOLD records
# Script by Rebecca E. Garner (rebecca.garner@mail.concordia.ca)

# Load libraries
library(tidyverse)

# Load palettes
source("scripts/00_palettes.R")

#### Import and format data ####
# Import BOLD records
bold_records <- read_tsv("output/2022-06-18_bold_records.tsv", col_names = TRUE, col_types = cols(.default = "c"))

# Extract taxonomy and conservation status information
conservation_status <- bold_records %>%
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


#### Examine conservation statuses of freshwater taxa ####
# Number of species per conservation status category
(n_species_conservation <- conservation_status %>%
   group_by(taxon, assessment) %>%
   count(name = "n_species") %>%
   ungroup() %>%
   pivot_wider(names_from = taxon, values_from = n_species, values_fill = 0) %>%
   arrange(factor(assessment, levels = cosewic_ranking)) %>%
   mutate(n_species_total = rowSums(across(where(is.numeric)))))

# Determine number of taxa by conservation status
(n_taxa_by_assessment_desforges <- conservation_status %>%
    group_by(assessment) %>%
    count(name = "n_taxa_desforges") %>%
    ungroup() %>%
    arrange(-n_taxa_desforges))


#### Examine BOLD records ####
# Determine number of total vs. unique BOLD records
bold_records %>%
  filter(!is.na(processid)) %>%
  nrow()  # Number of total BOLD records

bold_records %>%
  filter(!is.na(processid)) %>%
  distinct(processid, .keep_all = TRUE) %>%
  nrow()  # Number of unique BOLD records

# Determine species without records
(tax_without_records <- bold_records %>%
    filter(is.na(processid)) %>%
    distinct(species) %>%
    pull(species))
length(tax_without_records)  # Number of species without records

# Determine species with records
taxa_with_records <- bold_records %>%
  filter(!is.na(processid)) %>%
  distinct(species) %>%
  pull(species)
length(taxa_with_records)  # Number of species with records

# Count total BOLD records by species
n_bold_records <- bold_records %>%
  mutate(record = case_when(is.na(processid) ~ 0,
                            TRUE ~ 1)) %>%
  group_by(taxon, class, order, family, genus, species) %>%
  summarize(n_bold_records = sum(record),
            bold_processid = toString(processid)) %>%
  ungroup()

# Plot histogram of n BOLD records by taxon
(n_bold_records %>%
    left_join(conservation_status, by = "species") %>%
    ggplot() +
    geom_histogram(aes(x = n_bold_records, y = ..count..,
                       fill = factor(assessment, levels = cosewic_ranking)),
                   boundary = 0, binwidth = 10) +
    scale_fill_manual(values = palette_assessment) +
    labs(x = "Number of BOLD records",
         y = "Number of taxa",
         fill = "Assessment") +
    theme_bw())

(n_bold_records %>%
    left_join(conservation_status, by = "species") %>%
    ggplot() +
    facet_wrap(~taxon, scales = "free", nrow = 2) +
    geom_histogram(aes(x = n_bold_records, y = ..count..,
                       fill = factor(assessment, levels = cosewic_ranking)),
                   boundary = 0) +
    scale_fill_manual(values = palette_assessment) +
    labs(x = "Number of BOLD records",
         y = "Number of taxa",
         fill = "Assessment") +
    theme_bw())


#### Examine markers and primers in BOLD records ####
# Determine number of BOLD records with marker information
(nrecords_marker <- bold_records %>%
   filter(!is.na(processid)) %>%
   distinct(processid, .keep_all = TRUE) %>%
   filter(!is.na(markercode) | !is.na(marker_codes)) %>%
   nrow())

# Determine diversity of markers
markers <- bold_records %>%
  filter(!is.na(processid)) %>%
  distinct(processid, .keep_all = TRUE) %>%
  mutate(marker = case_when(!is.na(markercode) ~ markercode,
                            is.na(markercode) & !is.na(marker_codes) ~ marker_codes)) %>%
  group_by(marker) %>%
  count(name = "n_records") %>%
  ungroup()

# (Collapse markers with repetitive names)
markers_abbr <- markers %>%
  filter(!marker %in% unique(bold_records$markercode)) %>%
  mutate(marker_split = strsplit(as.character(marker), "\\|")) %>%
  unnest(marker_split) %>%
  distinct(marker, n_records, marker_split) %>%
  group_by(marker, n_records) %>%
  summarize(marker_abbr = toString(marker_split)) %>%
  ungroup() %>%
  mutate(marker_abbr = str_replace_all(marker_abbr, ", ", "|")) %>%
  group_by(marker_abbr) %>%
  summarize(n_records = sum(n_records)) %>%
  ungroup() %>%
  rename(marker = marker_abbr)

# (Join all markers counts)
markers %>%
  filter(marker %in% unique(bold_records$markercode)) %>%
  bind_rows(markers_abbr) %>%
  group_by(marker) %>%
  summarize(n_records = sum(n_records)) %>%
  ungroup() %>%
  arrange(-n_records) %>%
  filter(!is.na(marker)) %>%
  mutate(pct_records = n_records/nrecords_marker * 100)


#### Examine BOLD records with sequence information ####
# Calculate sequence length
seq_lengths <- bold_records %>%
  filter(!is.na(processid)) %>%
  distinct(processid, .keep_all = TRUE) %>%
  distinct(nucleotides) %>%
  filter(!is.na(nucleotides)) %>%
  mutate(dna_unaligned = str_replace_all(nucleotides, "-", "")) %>%
  mutate(seq_length = nchar(dna_unaligned))

(seqlengths_histogram <- seq_lengths %>%
    ggplot() +
    geom_histogram(aes(x = seq_length, y = ..count..), binwidth = 10) +
    labs(x = "Sequence length (bp)",
         y = "Number of barcodes") +
    theme_bw())

mean(seq_lengths$seq_length)
range(seq_lengths$seq_length)


#### Examine BOLD records by taxonomic group ####
# Determine number of taxa with BOLD records by "taxon"
(n_taxa_by_taxon_bold <- bold_records %>%
   filter(!is.na(processid)) %>%
   distinct(taxon, species) %>%
   group_by(taxon) %>%
   count(name = "n_taxa_bold") %>%
   ungroup() %>%
   arrange(-n_taxa_bold))
sum(n_taxa_by_taxon_bold$n_taxa_bold) == length(taxa_with_records)  # Should evaluate to TRUE

# Determine number of taxa with BOLD records by order
(n_taxa_by_order_bold <- bold_records %>%
    filter(!is.na(processid)) %>%
    distinct(taxon, class, order, species) %>%
    group_by(taxon, class, order) %>%
    count(name = "n_taxa_bold") %>%
    ungroup() %>%
    arrange(-n_taxa_bold))
sum(n_taxa_by_order_bold$n_taxa_bold) == length(taxa_with_records)  # Should evaluate to TRUE

# Summarize BOLD records by "taxon"
(bold_records_summary_by_taxon <- bold_records %>%
    filter(!is.na(processid)) %>%
    group_by(taxon) %>%
    count(name = "n_records") %>%
    ungroup() %>%
    arrange(-n_records) %>%
    left_join(n_taxa_by_taxon_bold, "taxon") %>%
    left_join(n_taxa_by_taxon_desforges, "taxon") %>%
    mutate(pct_taxa_bold = n_taxa_bold/n_taxa_desforges * 100))

# Summarize BOLD records by order
(bold_records_summary_by_order <- bold_records %>%
    filter(!is.na(processid)) %>%
    group_by(taxon, class, order) %>%
    count(name = "n_records") %>%
    ungroup() %>%
    arrange(-n_records) %>%
    left_join(n_taxa_by_order_bold, c("taxon", "class", "order")) %>%
    left_join(n_taxa_by_order_desforges, c("taxon", "class", "order")) %>%
    mutate(pct_taxa_bold = n_taxa_bold/n_taxa_desforges * 100))

# Summarize BOLD records by species
(bold_records_summary_by_species <- bold_records %>%
    filter(!is.na(processid)) %>%
    group_by(taxon, class, order, family, genus, species) %>%
    count(name = "n_records") %>%
    ungroup())


#### Examine BOLD records by conservation status ####
# Determine BOLD records by conservation status assessment
(bold_records_summary_by_assessment <- bold_records %>%
   filter(!is.na(processid)) %>%
   distinct(taxon, class, order, family, genus, species, assessment) %>%
   group_by(assessment) %>%
   count(name = "n_taxa") %>%
   ungroup() %>%
   left_join(n_taxa_by_assessment_desforges, "assessment") %>%
   mutate(pct_taxa_bold = n_taxa/n_taxa_desforges * 100) %>%
   arrange(-pct_taxa_bold))
