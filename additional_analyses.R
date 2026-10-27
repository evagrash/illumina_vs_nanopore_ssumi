################################################################################
#                                                                              #
#                        FOR REVIEW - ADDITIONAL ANALYSES                      #
#                                                                              #
################################################################################

# source("~/R/nanopore_vs_illumina_2025/analyses.R")
# Loads required objects

########################################################################################
# Reviewer 1 Major comment 4: Number of sequences (library size) per sample and method #
########################################################################################

# Calculated before filtering, before decontamiantion 

library_sizes_raw <- bind_rows(
  
  v3v4_otutab %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    mutate(Library_Size = rowSums(across(where(is.numeric))),
           Method = "V3-V4 Illumina") %>%
    select(Sample, Library_Size, Method),
  
  v1v3_otutab %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    mutate(Library_Size = rowSums(across(where(is.numeric))),
           Method = "V1-V3 Illumina") %>%
    select(Sample, Library_Size, Method),
  
  ssumi_otutab %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    mutate(Library_Size = rowSums(across(where(is.numeric))),
           Method = "ssUMI V1-V9 ONT") %>%
    select(Sample, Library_Size, Method)

) %>%
  filter(!Sample %in% exclude_samples) %>%
  left_join(metadata, by = "Sample") %>%
  filter(!is.na(sample_mock_control))

# Mock averages
library_sizes_raw %>%
  filter(sample_mock_control == "mock") %>%
  group_by(Method, Person) %>%
  summarise(
    mean_reads = mean(Library_Size),
    sd_reads   = sd(Library_Size),
    min_reads  = min(Library_Size),
    max_reads  = max(Library_Size),
    n          = n(),
    .groups    = "drop") %>%
  arrange(Method, Person) %>%
  print(n = 30)

# Sample averages
library_sizes_raw %>%
  filter(sample_mock_control == "sample") %>%
  group_by(Method) %>%
  summarise(
    mean_reads = mean(Library_Size),
    sd_reads   = sd(Library_Size),
    min_reads  = min(Library_Size),
    max_reads  = max(Library_Size),
    n          = n(),
    .groups    = "drop") %>%
  arrange(Method) %>%
  print()

################################################################################
# Reviewer 1 Major comment 4: UMI coverage per sample                          #
################################################################################

umi_counts %>%
  group_by(Sample_Name) %>%
  summarise(
    mean_coverage = mean(ubs_count),
    sd_coverage   = sd(ubs_count),
    median_coverage = median(ubs_count),
    n_umis        = n(),
    .groups = "drop") %>%
  arrange(Sample_Name) %>%
  print(n = 30)

umi_counts %>%
  filter(Sample_Name != "NK") %>%
  group_by(Sample_Name) %>%
  summarise(mean_per_sample = mean(ubs_count), .groups = "drop") %>%
  mutate(group = case_when(
    grepl("skin", Sample_Name, ignore.case = TRUE)  ~ "skin mock",
    grepl("staph", Sample_Name, ignore.case = TRUE) ~ "staph mock",
    TRUE ~ "sample")) %>%
  group_by(group) %>%
  summarise(
    mean_coverage = mean(mean_per_sample),
    sd_coverage   = sd(mean_per_sample),
    .groups = "drop")

################################################################################
# Reviewer 1 Major comment 4: Negative control genera                          #
################################################################################

physeq_v1v3_tax <- phyloseq(
  otu_table(as.matrix(v1v3_otutab), taxa_are_rows = TRUE),
  sample_data(metadata %>% column_to_rownames("Sample")),
  tax_table(v1v3_vsearch_tax %>% column_to_rownames("OTU") %>% as.matrix()))

physeq_v3v4_tax <- phyloseq(
  otu_table(as.matrix(v3v4_otutab), taxa_are_rows = TRUE),
  sample_data(metadata %>% column_to_rownames("Sample")),
  tax_table(v3v4_vsearch_tax %>% column_to_rownames("OTU") %>% as.matrix()))

physeq_ssumi_tax <- phyloseq(
  otu_table(as.matrix(ssumi_otutab), taxa_are_rows = TRUE),
  sample_data(metadata %>% column_to_rownames("Sample")),
  tax_table(ssumi_vsearch_tax %>% column_to_rownames("OTU") %>% as.matrix()))

# V3-V4 negative controls
prune_samples(c("NPCP_NKPa2", "NPCP_NKPa"), physeq_v3v4_tax) %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  group_by(Sample, Genus) %>%
  summarise(reads = sum(Abundance), .groups = "drop") %>%
  arrange(Sample, desc(reads)) %>%
  print(n = 50)

# V1-V3 negative controls
prune_samples(c("NPCP_NKPb2", "NPCP_NKPb"), physeq_v1v3_tax) %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  group_by(Sample, Genus) %>%
  summarise(reads = sum(Abundance), .groups = "drop") %>%
  arrange(Sample, desc(reads)) %>%
  print(n = 50)

# ssUMI negative control
prune_samples(c("NK"), physeq_ssumi_tax) %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  group_by(Sample, Genus) %>%
  summarise(reads = sum(Abundance), .groups = "drop") %>%
  arrange(Sample, desc(reads)) %>%
  print(n = 50)

################################################################################
# Reviewer 1 Major comment 2: L1 distance to theoretical, with uncertainity    #                                      
################################################################################

# SKIN L1 WITH ST.DEV ON SPECIES LEVEL -----------------------------------------

l1_plot_data_wide <- l1_plot_data_replicates %>%
  select(-Sample) %>%
  pivot_wider(names_from = Primers,
              values_from = rel_abund,
              values_fill = 0)

theoretical_rel <- theoretical_plus_melt %>%
  filter(Primers == "Theoretical") %>%
  group_by(Species2) %>%
  summarise(theoretical = sum(Abundance), .groups = "drop") %>%
  mutate(
    theoretical = theoretical / sum(theoretical),
    Species2 = gsub("_", " ", Species2))

l1_plot_data_wide %>%
  left_join(theoretical_rel, by = "Species2") %>%
  mutate(theoretical = replace_na(theoretical, 0)) %>%  # (other) categories get abundance == 0 in theoretical mock
  pivot_longer(cols = c(V1V3, V1V9, V3V4), names_to = "Method", values_to = "rel_abund") %>%
  group_by(Method, Replicate) %>%
  summarise(
    L1 = sum(abs(rel_abund - theoretical), na.rm = TRUE),
    .groups = "drop") %>%
  group_by(Method) %>%
  summarise(
    mean_L1 = mean(L1),
    sd_L1   = sd(L1),
    .groups = "drop") %>%
  mutate(label = sprintf("%.2f ± %.3f", mean_L1, sd_L1))

skin_l1_per_replicate <- l1_plot_data_wide %>%
  left_join(theoretical_rel, by = "Species2") %>%
  mutate(theoretical = replace_na(theoretical, 0)) %>%
  pivot_longer(cols = c(V1V3, V1V9, V3V4), names_to = "Method", values_to = "rel_abund") %>%
  group_by(Method, Replicate) %>%
  summarise(L1 = sum(abs(rel_abund - theoretical)), .groups = "drop")

# STAPH L1 WITH ST.DEV ON SPECIES LEVEL ----------------------------------------

staph_l1_plot_data_wide <- staph_l1_plot_data_replicates %>%
  select(-Sample) %>%
  pivot_wider(names_from = Primers,
              values_from = rel_abund,
              values_fill = 0)

staph_theoretical_rel <- staph_theoretical_plus_melt %>%
  filter(Primers == "Theoretical") %>%
  group_by(Species2) %>%
  summarise(theoretical = sum(Abundance), .groups = "drop") %>%
  mutate(
    theoretical = theoretical / sum(theoretical),
    Species2 = gsub("_", " ", Species2))

staph_l1_plot_data_wide %>%
  left_join(staph_theoretical_rel, by = "Species2") %>%
  mutate(theoretical = replace_na(theoretical, 0)) %>%
  pivot_longer(cols = c(V1V3, V1V9, V3V4), names_to = "Method", values_to = "rel_abund") %>%
  group_by(Method, Replicate) %>%
  summarise(
    L1 = sum(abs(rel_abund - theoretical), na.rm = TRUE),
    .groups = "drop") %>%
  group_by(Method) %>%
  summarise(
    mean_L1 = mean(L1),
    sd_L1   = sd(L1),
    .groups = "drop") %>%
  mutate(label = sprintf("%.2f ± %.3f", mean_L1, sd_L1))

staph_l1_per_replicate <- staph_l1_plot_data_wide %>%
  left_join(staph_theoretical_rel, by = "Species2") %>%
  mutate(theoretical = replace_na(theoretical, 0)) %>%
  pivot_longer(cols = c(V1V3, V1V9, V3V4), names_to = "Method", values_to = "rel_abund") %>%
  group_by(Method, Replicate) %>%
  summarise(L1 = sum(abs(rel_abund - theoretical)), .groups = "drop")

# SKIN L1 WITH STDEV ON GENUS LEVEL AND PLOT -----------------------------------

# theoretical at genus level
theoretical_genus_rel <- theoretical_plus_melt %>%
  filter(Primers == "Theoretical") %>%
  mutate(Genus2 = gsub("_", " ", Genus2)) %>%
  group_by(Genus2) %>%
  summarise(theoretical = sum(Abundance), .groups = "drop") %>%
  mutate(theoretical = theoretical / sum(theoretical))

# observed averaged across replicates, joined with theoretical
genus_vs_truth <- l1_genus_wide %>%
  group_by(Genus2) %>%
  summarise(V1V9 = mean(V1V9), V1V3 = mean(V1V3), V3V4 = mean(V3V4), .groups = "drop") %>%
  left_join(theoretical_genus_rel, by = "Genus2") %>%
  mutate(theoretical = replace_na(theoretical, 0)) %>%
  pivot_longer(cols = c(V1V3, V1V9, V3V4), names_to = "Method", values_to = "observed") %>%
  mutate(Method = factor(Method, levels = c("V1V9", "V3V4", "V1V3")))

# skin L1 on genus level calculation
skin_l1_per_replicate_genus <- l1_genus_wide %>%
  left_join(theoretical_genus_rel, by = "Genus2") %>%
  mutate(theoretical = replace_na(theoretical, 0)) %>%
  pivot_longer(cols = c(V1V3, V1V9, V3V4), names_to = "Method", values_to = "rel_abund") %>%
  group_by(Method, Replicate) %>%
  summarise(L1 = sum(abs(rel_abund - theoretical)), .groups = "drop")

skin_l1_per_replicate_genus %>%
  group_by(Method) %>%
  summarise(
    mean_L1 = mean(L1),
    sd_L1   = sd(L1),
    .groups = "drop")

# PLOT WITH ERROR BARS ---------------------------------------------------------

# difference between observed and theoretical
genus_deviation <- genus_vs_truth %>%
  mutate(deviation = observed - theoretical)

# three theoretical abundance values for vertical lines
theoretical_vals <- genus_deviation %>%
  distinct(Genus2, theoretical) %>%
  pull(theoretical) %>%
  unique() %>%
  sort()

# per-replicate genus deviation
genus_deviation_reps <- l1_genus_wide %>%
  left_join(theoretical_genus_rel, by = "Genus2") %>%
  mutate(theoretical = replace_na(theoretical, 0)) %>%
  pivot_longer(cols = c(V1V3, V1V9, V3V4), names_to = "Method", values_to = "observed") %>%
  mutate(
    deviation = observed - theoretical,
    Method = factor(Method, levels = c("V1V9", "V3V4", "V1V3")))

# mean and SD of deviation across replicates
genus_deviation_avg <- genus_deviation_reps %>%
  group_by(Genus2, Method, theoretical) %>%
  summarise(
    mean_dev = mean(deviation),
    sd_dev   = sd(deviation),
    .groups  = "drop")

# order alphabetically, A on top
genus_deviation_avg <- genus_deviation_avg %>%
  mutate(Genus2 = factor(Genus2, levels = rev(sort(unique(as.character(Genus2))))))

genus_deviation_avg <- genus_deviation_avg %>%
  mutate(Genus2 = factor(Genus2, levels = rev(genus_order)))

# plot with error bars
new_l1_genus_to_truth <- ggplot(genus_deviation_avg) +
  geom_vline(xintercept = theoretical_vals,
             color = "grey70", linewidth = 0.4, linetype = "dashed") +
  geom_segment(aes(x = theoretical,
                   xend = theoretical + mean_dev,
                   y = Genus2,
                   yend = Genus2),
               color = "black", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = theoretical + mean_dev - sd_dev,
                     xmax = theoretical + mean_dev + sd_dev,
                     y = Genus2),
                 height = 0.5, linewidth = 0.4) +
  geom_point(aes(x = theoretical + mean_dev,
                 y = Genus2,
                 fill = Genus2),
             shape = 21, size = 2.5, color = "black") +
  facet_wrap(~ Method, ncol = 3,
             labeller = labeller(Method = c(
               "V1V9" = "ssUMI V1-V9 ONT",
               "V3V4" = "V3-V4 Illumina",
               "V1V3" = "V1-V3 Illumina"))) +
  scale_fill_manual(values = genus_colors, name = "Genus") +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 0.5)) +
  labs(x = "Relative abundance", y = NULL) +
  theme_bw() +
  theme(
    axis.text.y = element_text(face = "italic"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none")

ggsave("~/R/nanopore_vs_illumina_2025/Paper_images_ideas/new_l1_genus_to_truth.svg", new_l1_genus_to_truth, width = 170, height = 90, units = "mm", dpi = 600)
ggsave("~/R/nanopore_vs_illumina_2025/Paper_images_ideas/new_l1_genus_to_truth.png", new_l1_genus_to_truth, width = 170, height = 90, units = "mm", dpi = 600)

################################################################################
# Reviewer 2: Minor comments (104-106, 104): Recalculate how many reads were   # 
# classified at a species level in eack mock x primer pair                     #
################################################################################

how_many_other_reps <- theoretical_plus_melt %>%
  filter(Primers != "Theoretical",
         sample_mock_control == "mock") %>%
  group_by(Primers, Sample, Species2) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Primers, Sample) %>%
  mutate(Percentage = Abundance / sum(Abundance) * 100) %>%
  ungroup() %>%
  mutate(is_other = str_detect(Species2, "\\(other\\)")) %>% 
  group_by(Primers, Sample) %>%
  summarise(
    pct_species_classified = sum(Percentage[!is_other]),
    .groups = "drop") %>%
  group_by(Primers) %>%
  summarise(
    mean_pct = mean(pct_species_classified),
    sd_pct   = sd(pct_species_classified),
    .groups  = "drop") %>%
  mutate(label = sprintf("%.1f ± %.1f%%", mean_pct, sd_pct))

how_many_other_reps_staph <- staph_theoretical_plus_melt %>%
  filter(Primers != "Theoretical",
         sample_mock_control == "mock") %>%
  group_by(Primers, Sample, Species2) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Primers, Sample) %>%
  mutate(Percentage = Abundance / sum(Abundance) * 100) %>%
  ungroup() %>%
  mutate(is_other = grepl("\\(other\\)", Species2)) %>%
  group_by(Primers, Sample) %>%
  summarise(
    pct_species_classified = sum(Percentage[!is_other]),
    .groups = "drop") %>%
  group_by(Primers) %>%
  summarise(
    mean_pct = mean(pct_species_classified),
    sd_pct   = sd(pct_species_classified),
    .groups  = "drop") %>%
  mutate(label = sprintf("%.1f ± %.1f%%", mean_pct, sd_pct))

####################################################################################
# Reviewer 2 Minor comments: table of statistics from sequencing, raw ZOTU counts  #
####################################################################################

mean_zotus <- function(physeq, person_val) {
  keep <- sample_data(physeq)$Person == person_val
  ps <- prune_samples(keep, physeq)
  otu_mat <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)
  per_sample <- colSums(otu_mat > 0)
  list(mean = round(mean(per_sample), 1), sd = round(sd(per_sample), 1))
}

zotu_count_skin <- tibble(
  Primers = c("V1V9", "V3V4", "V1V3"),
  mean_zotus = c(
    mean_zotus(physeq_ssumi_tax, "skin_mock")$mean,
    mean_zotus(physeq_v3v4_tax,  "skin_mock")$mean,
    mean_zotus(physeq_v1v3_tax,  "skin_mock")$mean),
  sd_zotus = c(
    mean_zotus(physeq_ssumi_tax, "skin_mock")$sd,
    mean_zotus(physeq_v3v4_tax,  "skin_mock")$sd,
    mean_zotus(physeq_v1v3_tax,  "skin_mock")$sd))

zotu_count_staph <- tibble(
  Primers = c("V1V9", "V3V4", "V1V3"),
  mean_zotus = c(
    mean_zotus(physeq_ssumi_tax, "staph_mock")$mean,
    mean_zotus(physeq_v3v4_tax,  "staph_mock")$mean,
    mean_zotus(physeq_v1v3_tax,  "staph_mock")$mean),
  sd_zotus = c(
    mean_zotus(physeq_ssumi_tax, "staph_mock")$sd,
    mean_zotus(physeq_v3v4_tax,  "staph_mock")$sd,
    mean_zotus(physeq_v1v3_tax,  "staph_mock")$sd))

################################################################################
# Both reviewers: Missing species.                                             #
################################################################################

# SKIN MOCK - species in theoretical but missing from each method
skin_theoretical_species <- theoretical_plus_melt %>%
  filter(Primers == "Theoretical", Abundance > 0) %>%
  mutate(Species2 = gsub("_", " ", Species2)) %>%
  filter(Species2 != "Other") %>%
  pull(Species2) %>%
  unique()

skin_detected <- theoretical_plus_melt %>%
  filter(Primers != "Theoretical",
         sample_mock_control == "mock") %>%
  mutate(Species2 = gsub("_", " ", Species2)) %>%
  filter(!grepl("Other", Species2)) %>%
  group_by(Primers, Species2) %>%
  summarise(total = sum(Abundance), .groups = "drop") %>%
  filter(total > 0)

skin_missing <- tibble(Primers = c("V1V9", "V3V4", "V1V3")) %>%
  mutate(missing = map(Primers, ~{
    detected <- skin_detected %>% filter(Primers == .x) %>% pull(Species2)
    setdiff(skin_theoretical_species, detected)
  }))

# STAPH MOCK - species in theoretical but missing from each method
staph_theoretical_species <- staph_theoretical_plus_melt %>%
  filter(Primers == "Theoretical", Abundance > 0) %>%
  filter(!grepl("\\(other\\)", Species2)) %>%
  pull(Species2) %>%
  unique()

staph_detected <- staph_theoretical_plus_melt %>%
  filter(Primers != "Theoretical",
         sample_mock_control == "mock") %>%
  filter(!grepl("\\(other\\)", Species2)) %>%
  group_by(Primers, Species2) %>%
  summarise(total = sum(Abundance), .groups = "drop") %>%
  filter(total > 0)

staph_missing <- tibble(Primers = c("V1V9", "V3V4", "V1V3")) %>%
  mutate(missing = map(Primers, ~{
    detected <- staph_detected %>% filter(Primers == .x) %>% pull(Species2)
    setdiff(staph_theoretical_species, detected)
  }))

################################################################################
# Reviewer 1: Major comment 4: Add centroid positions to plot                  #
################################################################################

samples_beta_dispersion_n <- ggplot(all_bd_distances, 
                                  aes(x = factor(Primers, levels = c("V1V9", "V3V4", "V1V3")), 
                                      y = Dispersion)) +
  geom_line(aes(group = Person), color = "grey", linetype = "dashed") +
  geom_point(size = 4, aes(shape = Person)) +
  scale_shape_manual(values = c("A" = 2, "B" = 1, "C" = 4, "D" = 3, "E" = 5)) +
  theme_bw() +
  labs(y = "Distance to centroid on ZOTU level", x = "", title = "") +
  scale_x_discrete(labels = c("V1V9" = "ssUMI\nV1-V9\nONT",
                              "V3V4" = "V3-V4\nIllumina",
                              "V1V3" = "V1-V3\nIllumina")) +
  coord_cartesian(ylim = c(0, 0.6)) +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 10, margin = margin(t = 10)),
        axis.title.y = element_text(size = 10, margin = margin(r = 10)),
        axis.text  = element_text(size = 8),
        panel.border = element_rect(linewidth = 0.5))

genus_samples_beta_dispersion_n <- ggplot(genus_all_bd_distances, 
                                        aes(x = factor(Primers, levels = c("V1V9", "V3V4", "V1V3")), 
                                            y = Dispersion)) +
  geom_line(aes(group = Person), color = "grey", linetype = "dashed") +
  geom_point(size = 4, aes(shape = Person)) +
  scale_shape_manual(values = c("A" = 2, "B" = 1, "C" = 4, "D" = 3, "E" = 5)) +
  theme_bw() +
  labs(y = "Distance to centroid on genus level", x = "", title = "") +
  scale_x_discrete(labels = c("V1V9" = "ssUMI\nV1-V9\nONT",
                              "V3V4" = "V3-V4\nIllumina",
                              "V1V3" = "V1-V3\nIllumina")) +
  coord_cartesian(ylim = c(0, 0.6)) +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 10, margin = margin(t = 10)),
        axis.title.y = element_text(size = 10, margin = margin(r = 10)),
        axis.text  = element_text(size = 8),
        panel.border = element_rect(linewidth = 0.5))


combined_plot_n <- cowplot::plot_grid(samples_beta_dispersion_n, genus_samples_beta_dispersion_n, ncol = 2, labels = c("A", "B"))

ggsave("~/R/nanopore_vs_illumina_2025/Paper_images_ideas/beta_dispersion_combined_new.pdf", combined_plot_n, width = 170, height = 110, units = "mm", dpi = 600)
ggsave("~/R/nanopore_vs_illumina_2025/Paper_images_ideas/beta_dispersion_combined_new.svg", combined_plot_n, width = 170, height = 110, units = "mm", dpi = 600)

################################################################################
# Reviewer 2: Major comment 1: How many other, what are most common ZOTUs      #
################################################################################

other_count <- theoretical_plus_melt %>%
  filter(str_detect(Species2, "\\(other\\)")) %>% 
  filter(Abundance > 0)

other_summary <- other_count %>%
  group_by(Primers, OTU, Genus2) %>%
  summarise(
    mean_abundance = mean(Abundance),
    sd_abundance = sd(Abundance),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  arrange(Primers, desc(mean_abundance))

#################################################################################
# Reviewer 2: Major comment 2: Staphylococcus distribution in real life samples #
#################################################################################

# RE-MELT DATA  -----------------------------------------------------------------

melt_v3v4_bar <- psmelt(rarefied_physeq_v3v4) %>%
  group_by(Sample, Person, Species, Genus, Phylum, Primers) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(rel_abund = Abundance / sum(Abundance)) %>%
  ungroup()

melt_v1v3_bar <- psmelt(rarefied_physeq_v1v3) %>%
  group_by(Sample, Person, Species, Genus, Phylum, Primers) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(rel_abund = Abundance / sum(Abundance)) %>%
  ungroup()

melt_ssumi_bar <- psmelt(rarefied_physeq_ssumi) %>%
  group_by(Sample, Person, Species, Genus, Phylum, Primers) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(rel_abund = Abundance / sum(Abundance)) %>%
  ungroup()

# SSTAPH SPECIES THAT WERE DETECTED? -------------------------------------------

staph_species <- bind_rows(
  melt_v3v4_bar,
  melt_v1v3_bar,
  melt_ssumi_bar) %>%
  filter(Genus == "Staphylococcus", !is.na(Species)) %>%
  filter(Person %in% c("A", "B", "C", "D", "E")) %>%
  group_by(Primers, Species) %>%
  summarise(total = sum(Abundance), .groups = "drop") %>%
  group_by(Primers) %>%
  ungroup()

# Check detected to write a list of most common species
print(staph_species, n = 50)

# BIND AND ASSIGN DISPLAY CATEGORY --------------------------------------------

skin_phyla <- c("Firmicutes", "Actinobacteriota", "Proteobacteria", "Bacteroidota")

staph_species_show <- c(
  "Staphylococcus_epidermidis",
  "Staphylococcus_hominis",
  "Staphylococcus_capitis")

samples_bar_long <- bind_rows(melt_v3v4_bar, melt_v1v3_bar, melt_ssumi_bar) %>%
  filter(Person %in% c("G", "L", "M", "N", "U")) %>%
  mutate(
    display_cat = case_when(
        Genus == "Staphylococcus" & Species %in% staph_species_show ~ Species,
        Genus == "Staphylococcus" & Species %in% c("Staphylococcus_sp.", "uncultured_bacterium") ~ "Unidentified Staphylococcus",
        Genus == "Staphylococcus" ~ "Other Staphylococcus",
        Genus == "Cutibacterium" ~ "Cutibacterium",
        Phylum == "Firmicutes" ~ "Other Firmicutes",
        Phylum == "Actinobacteriota" ~ "Other Actinobacteriota",
        Phylum %in% skin_phyla ~ Phylum,
        TRUE ~ "Other"))

# AGGREGATE -------------------------------------------------------------------

bar_data <- samples_bar_long %>%
  group_by(Person, Primers, Sample, display_cat) %>%
  summarise(rel_abund = sum(rel_abund), .groups = "drop")

# FACTOR ORDER
cat_order <- c(
  "Staphylococcus_epidermidis",
  "Staphylococcus_capitis",
  "Staphylococcus_hominis",
  "Other Staphylococcus",
  "Unidentified Staphylococcus",
  "Other Firmicutes",
  "Cutibacterium",
  "Other Actinobacteriota",
  "Proteobacteria",
  "Bacteroidota",
  "Deinococcota",
  "Other")

bar_data <- bar_data %>%
  mutate(
    display_cat = factor(display_cat, levels = cat_order),
    Primers = factor(Primers, levels = c("V1V9", "V3V4", "V1V3"),
                     labels = c("ssUMI V1-V9 ONT", "V3-V4 Illumina", "V1-V3 Illumina")))

# COLORS ----------------------------------------------------------------------

bar_colors <- c(
  # Staphylococcus
  "Staphylococcus_auricularis"   = "#FDD49E",
  "Staphylococcus_capitis"       = "#FFE873",
  "Staphylococcus_epidermidis"   = "#FE9929",
  "Staphylococcus_haemolyticus"  = "#FDAE6B",
  "Staphylococcus_hominis"       = "#FEC44F",
  "Staphylococcus_lugdunensis"   = "#FFF0C0",
  "Staphylococcus_pasteuri"      = "#FFE873",
  "Unidentified Staphylococcus"  = "#7F2704",
  "Other Staphylococcus"         = "#E8601A", 
  
  # Firmicutes 
  "Other Firmicutes"             = "#CB181D",
  
  # Cutibacterium 
  "Cutibacterium"                = "#054A7A",
  
  # Actinomycetota 
  "Other Actinobacteriota"      = "#9ECAE1",
  
  # Other 
  "Proteobacteria"               = "#74C476",
  "Bacteroidota"                 = "#D9F0A3",
  "Deinococcota"                 = "#c9645d",
  "Other"                        = "#D9D9D9")

# PLOT ------------------------------------------------------------------------

phylum_bar <- ggplot(bar_data,
                     aes(x = Primers, y = rel_abund, fill = display_cat)) +
  geom_bar(stat = "identity", width = 0.8) +
  facet_grid(cols = vars(Person)) +
  scale_fill_manual(
    values = bar_colors,
    name   = NULL,
    guide  = guide_legend(reverse = FALSE)) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = c(0, 0)) +
  labs(x = NULL, y = "Relative abundance") +
  theme_bw() +
  theme(
    legend.position = "bottom", 
    legend.key.size = unit(0.25, "cm"), 
    legend.key.spacing.x = unit(0, "pt"),
    legend.key.spacing.y = unit(0, "pt"),
    legend.title.position = "top",
    legend.box = "vertical",
    axis.title = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text  = element_text(size = 8),
    legend.text = element_text(size = 6.3),
    legend.title = element_text(size = 9),
    panel.border = element_rect(linewidth = 0.5)) +
  guides(fill = guide_legend(byrow = FALSE, ncol = 2))

ggsave("~/R/nanopore_vs_illumina_2025/Paper_images_ideas/staphylococci_in_samples2.pdf", phylum_bar, width = 170, height = 150, units = "mm", dpi = 600)
ggsave("~/R/nanopore_vs_illumina_2025/Paper_images_ideas/staphylococci_in_samples2.svg", phylum_bar, width = 170, height = 150, units = "mm", dpi = 600)

#################################################################################################
# Reviewer 2, Major comment 3: Deseq2 (for staphylococci in skin mock community at genus level) #
#################################################################################################

library(DESeq2)

physeq_mock <- subset_samples(all_primers_phyloseq, Person == "skin_mock")

physeq_mock_genus <- tax_glom(physeq_mock, taxrank = "Genus", NArm = FALSE)

dds_mock <- phyloseq_to_deseq2(physeq_mock_genus, ~ Primers)
dds_mock <- DESeq(dds_mock, test = "Wald", fitType = "mean")

resultsNames(dds_mock) 

# Pairwise comparisons

contrasts_mock <- list(
  c("Primers", "V1V9", "V1V3"),
  c("Primers", "V3V4", "V1V3"),
  c("Primers", "V1V9", "V3V4"))

deseq_mock_results <- lapply(contrasts_mock, function(contrast) {
  res <- results(dds_mock, contrast = contrast, alpha = 0.05)
  as.data.frame(res) %>%
    rownames_to_column("OTU") %>%
    left_join(
      as.data.frame(tax_table(physeq_mock_genus)) %>% rownames_to_column("OTU"),
      by = "OTU") %>%
    mutate(comparison = paste(contrast[2], "vs", contrast[3])) %>%
    filter(!is.na(padj))}) %>%
  bind_rows()

deseq_mock_sig <- deseq_mock_results %>%
  filter(padj < 0.05) %>%
  select(comparison, Genus, log2FoldChange, padj) %>%
  arrange(comparison, padj)

deseq_mock_sig %>%
  filter(Genus %in% c("Staphylococcus", "Cutibacterium")) %>%
  select(comparison, Genus, log2FoldChange, padj) %>%
  arrange(Genus, comparison)

###################################################################################
# Reviewer 2 Minor comments: Shannon per samples (instead or ZOTU count/richness) #
###################################################################################

alpha_div %>%
  filter(sample_mock_control == "sample") %>%
  group_by(Primers) %>%
  summarise(
    mean_shannon = mean(Shannon),
    sd_shannon   = sd(Shannon),
    n            = n())


