################################################################################
# 1. FILTER GENOMES FOR INSILICO MOCK                                          #
################################################################################

library(tidyverse)

meta <- read_tsv("mock_community_genomes_metadata.tsv")

# Extract two-word species name (handles subspecies entries)
meta <- meta %>%
  mutate(species = str_extract(`Organism Name`, "^\\S+\\s+\\S+"))

keep_all_species <- c("Kocuria marina")

step1 <- meta %>%
  filter(
    species %in% keep_all_species | `ANI Check status` == "OK")

species_counts <- step1 %>%
  count(species)

print(species_counts, n = Inf)

# For species with > max_genomes, apply 4 contig filter
step2 <- step1 %>%
  group_by(species) %>%
  mutate(n_step1 = n()) %>%
  ungroup() %>%
  filter(
    n_step1 <= 30 |
      (n_step1 > 30 & `Assembly Stats Number of Contigs` %in% c(1, 2, 3, 4))
  )

species_counts <- step2 %>%
  count(species)

print(species_counts, n = Inf)

# Keep all if n <= max_genomes, randomly select max_genomes if n > max_genomes
step3 <- step2 %>%
  group_by(species) %>%
  group_modify(~ slice_sample(.x, n = min(nrow(.x), 40))) %>%
  ungroup()

species_counts <- step3 %>%
  count(species)

print(species_counts, n = Inf)

# Export
write_tsv(step3 %>% select("Assembly Accession"), "selected_accessions.tsv", col_names = FALSE)
write_tsv(step3, "selected_genomes_metadata.tsv")

###################################################################################
# 2. THE BASH PART - DOWNLOAD GENOMES, EXTRACT 16s, IN SILICO PCR, USEARCH, QIIME #
###################################################################################

################################################################################
# 3. IN SILICO MOCK COMMUNITY ANALYSIS - cleaned up                            #
################################################################################

library(tidyverse)
library(qiime2R)

# IMPORT -----------------------------------------------------------------------

# source("~/R/nanopore_vs_illumina_2025/analyses.R")
# source("~/R/nanopore_vs_illumina_2025/additional_analyses.R")

insilico_v1v3_otutab <- read_tsv("~/R/nanopore_vs_illumina_2025/in_silico_mock/V1V3_zotutab.tsv") %>% 
  column_to_rownames("#OTU ID")

insilico_v3v4_otutab <- read_tsv("~/R/nanopore_vs_illumina_2025/in_silico_mock/V3V4_zotutab.tsv") %>% 
  column_to_rownames("#OTU ID")

insilico_v1v9_otutab <- read_tsv("~/R/nanopore_vs_illumina_2025/in_silico_mock/V1V9_zotutab.tsv") %>% 
  column_to_rownames("#OTU ID")

import_qiime_taxonomy <- function(qza_taxonomy_path) {
  read_qza(qza_taxonomy_path)$data %>%
    rename(OTU = `Feature.ID`) %>%
    mutate(Taxon = str_replace_all(Taxon, "\\\\|\\\"|\\(\\d+\\)", ""),
           Taxon = str_replace(Taxon, ";$", "")) %>%
    separate(Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
             sep = ";", fill = "right") %>%
    mutate(across(Domain:Species, ~ str_trim(.)),
           across(Domain:Species, ~ str_remove(., "^[a-z]__")),
           across(Domain:Species, ~ str_replace_all(., "_", " ")),
           across(Domain:Species, ~ na_if(str_trim(.), "")))
}

insilico_v1v3_vtax <- import_qiime_taxonomy("~/R/nanopore_vs_illumina_2025/in_silico_mock/V1V3_vsearchZOTUclassification.qza")
insilico_v3v4_vtax <- import_qiime_taxonomy("~/R/nanopore_vs_illumina_2025/in_silico_mock/V3V4_vsearchZOTUclassification.qza")
insilico_v1v9_vtax <- import_qiime_taxonomy("~/R/nanopore_vs_illumina_2025/in_silico_mock/V1V9_vsearchZOTUclassification.qza")

# JOIN ZOTUS + TAXONOMY --------------------------------------------------------

read_zotu <- function(df, region) {
  df %>%
    rownames_to_column("ZOTU") %>%
    pivot_longer(-ZOTU, names_to = "Accession", values_to = "count") %>%
    filter(count > 0) %>%
    mutate(Primers = region)
}

tax_all <- bind_rows(
  insilico_v1v3_vtax %>% mutate(Primers = "V1V3"),
  insilico_v3v4_vtax %>% mutate(Primers = "V3V4"),
  insilico_v1v9_vtax %>% mutate(Primers = "V1V9"))

all_zotus <- bind_rows(
  read_zotu(insilico_v1v3_otutab, "V1V3"),
  read_zotu(insilico_v3v4_otutab, "V3V4"),
  read_zotu(insilico_v1v9_otutab, "V1V9")) %>%
  left_join(
    meta %>% mutate(Accession = str_remove(`Assembly Accession`, "\\.\\d+$")),
    by = "Accession") %>%
  left_join(tax_all, by = c("ZOTU" = "OTU", "Primers"))

# DETECTION TABLE IN SEQUENCED MOCKS -------------------------------------------

all_primers <- c("V1V3", "V1V9", "V3V4")

all_mock_data <- bind_rows(theoretical_plus_melt, staph_theoretical_plus_melt)

theoretical_species <- all_mock_data %>%
  filter(Sample == "Theoretical") %>%
  distinct(Person, Species)

detected_species <- all_mock_data %>%
  filter(Sample != "Theoretical", sample_mock_control == "mock") %>%
  group_by(Person, Primers, Species) %>%
  summarise(total_abundance = sum(Abundance), .groups = "drop") %>%
  filter(total_abundance > 0)

missed_species <- theoretical_species %>%
  crossing(Primers = unique(detected_species$Primers)) %>%
  anti_join(detected_species, by = c("Person", "Species", "Primers")) %>%
  filter(Primers != "Theoretical")
            
theoretical_combos <- bind_rows(
  theoretical_mock %>% mutate(species = str_replace_all(Species, "_", " ")),
  theoretical_staph_mock %>% mutate(species = str_replace_all(Species, "_", " "))) %>%
  distinct(species) %>%
  crossing(Primers = all_primers)

missed_any <- missed_species %>%
  mutate(species = str_replace_all(Species, "_", " ")) %>%
  group_by(species, Primers) %>%
  summarise(missed_real = TRUE, .groups = "drop")

detected_labels <- all_zotus %>%
  filter(!is.na(Species), !str_detect(Species, "uncultured|metagenome")) %>%
  transmute(Primers, species = Species) %>%
  distinct() %>%
  mutate(label_found = TRUE)

detection_table <- theoretical_combos %>%
  left_join(detected_labels, by = c("species", "Primers")) %>%
  left_join(missed_any, by = c("species", "Primers")) %>%
  mutate(
    insilico_detected = !is.na(label_found),
    detected_real = is.na(missed_real)) %>%
  select(species, Primers, insilico_detected, detected_real) %>%
  arrange(species, Primers)

detection_wide <- detection_table %>%
  pivot_wider(names_from = Primers, values_from = c(insilico_detected, detected_real)) %>%
  rename(
    insilico_V1V9 = insilico_detected_V1V9,
    insilico_V3V4 = insilico_detected_V3V4,
    insilico_V1V3 = insilico_detected_V1V3,
    mock_V1V9     = detected_real_V1V9,
    mock_V3V4     = detected_real_V3V4,
    mock_V1V3     = detected_real_V1V3)

# FINAL TABLE WITH READS AND DETECTION ########################################

# CLASSIFY ZOTUS DIRECTLY FROM all_zotus --------------------------------------

zotus_classified <- all_zotus %>%
  mutate(
    tax_cat = case_when(
      is.na(Species) | str_detect(Species, "uncultured|metagenome") ~ "genus",
      Species == species ~ "correct",
      TRUE ~ "wrong")) %>%
  group_by(species, Primers, tax_cat) %>%
  summarise(
    n_zotus = n_distinct(ZOTU),
    n_reads = sum(count),
    .groups = "drop") %>%
  pivot_wider(
    names_from  = tax_cat,
    values_from = c(n_zotus, n_reads),
    values_fill = 0) %>%
  mutate(
    across(any_of(c("n_zotus_correct", "n_zotus_wrong", "n_zotus_genus",
                    "n_reads_correct", "n_reads_wrong", "n_reads_genus")),
           ~ replace_na(., 0)),
    col_correct = paste0(n_zotus_correct, " (", n_reads_correct, ")"),
    col_wrong = paste0(n_zotus_wrong, " (", n_reads_wrong, ")"),
    col_genus = paste0(n_zotus_genus, " (", n_reads_genus, ")"))

# WIDE FORMAT ------------------------------------------------------------------

zotus_wide <- zotus_classified %>%
  select(species, Primers, col_correct, col_wrong, col_genus) %>%
  pivot_wider(
    names_from = Primers,
    values_from = c(col_correct, col_wrong, col_genus),
    values_fill = "0 (0)")

# FINAL TABLE ------------------------------------------------------------------

n_strains_true <- all_zotus %>%
  group_by(species) %>%
  summarise(No_strains = n_distinct(Accession), .groups = "drop")

final_table <- detection_wide %>%
  left_join(n_strains_true, by = "species") %>%
  left_join(zotus_wide, by = "species") %>%
  select(species, No_strains,
    col_correct_V1V9, col_wrong_V1V9, col_genus_V1V9, mock_V1V9,
    col_correct_V3V4, col_wrong_V3V4, col_genus_V3V4, mock_V3V4,
    col_correct_V1V3, col_wrong_V1V3, col_genus_V1V3, mock_V1V3) %>%
  mutate(across(starts_with("mock_"), ~ ifelse(., "Y", "N"))) %>%
  arrange(species)

write_csv2(final_table, "final_in_silico_mock_table.csv")

# SHARED ZOTUS -----------------------------------------------------------------  

shared_all <- all_zotus %>%
  group_by(Primers, ZOTU) %>%
  filter(n_distinct(species) > 1) %>%
  group_by(Primers, ZOTU, Species) %>%
  summarise(
    source_species = paste(unique(species), collapse = "; "),
    n_reads = sum(count),
    .groups = "drop") %>%
  group_by(Primers, ZOTU) %>%
  mutate(
    all_sources = paste(unique(source_species), collapse = "; "),
    n_species = n_distinct(source_species)) %>%
  ungroup() %>%
  select(Primers, ZOTU, all_sources, classified_as = Species, n_reads) %>%
  distinct() %>%
  arrange(Primers, ZOTU)

# ZOTUS PER SPECIES PER PRIMER WITH CLASSIFICATIONS ---------------------------

zotus_per_species <- all_zotus %>%
  mutate(
    display_tax = case_when(
      is.na(Species) ~ paste0(Genus, " NA"),
      str_detect(Species, "uncultured|metagenome") ~ paste0(Genus, " uncultured bacterium"),
      TRUE  ~ Species)) %>%
  group_by(species, Primers, ZOTU, display_tax) %>%
  summarise(n_strains_per_zotu = n_distinct(Accession), .groups = "drop") %>%
  group_by(species, Primers, display_tax) %>%
  summarise(
    n_zotus_this_tax = n_distinct(ZOTU),
    n_strains_this_tax = sum(n_strains_per_zotu),
    .groups = "drop") %>%
  group_by(species, Primers) %>%
  summarise(
    n_zotus = sum(n_zotus_this_tax),
    n_strains = sum(n_strains_this_tax),
    classifications = paste0(
      display_tax, " (", n_zotus_this_tax, " ZOTUs, ", n_strains_this_tax, "x)",
      collapse = "; "),
    .groups = "drop") %>%
  arrange(Primers, desc(n_zotus))

# MISCLASSIFIED SPECIES --------------------------------------------------------

misclassified <- zotus_per_species %>%
  rowwise() %>%
  mutate(
    wrong_species = {
      matches <- str_extract_all(classifications, "[A-Z][a-z]+ [a-z]+")[[1]]
      wrong <- matches[matches != species &
                       !str_detect(matches, "uncultured")]
      if (length(wrong) == 0) NA_character_ else paste(wrong, collapse = "; ")}) %>%
  ungroup() %>%
  filter(!is.na(wrong_species)) %>%
  select(species, Primers, wrong_species, classifications)

