################################################################################
#                                                                              #
#               ssUMI Nanopore V1-V9 vs Illumina V3-V4 and V1-V3               #
#                                                                              #
################################################################################

library(tidyverse)
library(vegan)
library(phyloseq)
library(ggrain)
library(ggpubr)
library(glue)
library(pairwiseAdonis)
library(emmeans)

renv::snapshot()

################################# IMPORT #######################################

# ZOTU tables

v1v3_otutab <- read_tsv("~/R/nanopore_vs_illumina_2025/input/zotutab_rawv1v3.tsv") %>% 
  column_to_rownames("#OTU ID")

v3v4_otutab <- read_tsv("~/R/nanopore_vs_illumina_2025/input/zotutab_rawv3v4.tsv") %>% 
  column_to_rownames("#OTU ID")

ssumi_otutab <- read_tsv("~/R/nanopore_vs_illumina_2025/input/ssumi_zotutab_raw.tsv") %>% 
  column_to_rownames("#OTU ID") %>%
  select(!c("barcode01", "unclassified_1", "unclassified_2"))

# Taxonomy tables

library(qiime2R)

import_qiime_taxonomy <- function(qza_taxonomy_path) {
  read_qza(qza_taxonomy_path)$data %>%
    rename(OTU = `Feature.ID`) %>%
    mutate(Taxon = str_replace_all(Taxon, "\\\\|\\\"|\\(\\d+\\)", ""),
           Taxon = str_replace(Taxon, ";$", "")) %>%
    separate(Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", fill = "right") %>%
    mutate(across(everything(), ~ gsub("^ [a-z]__*", "", .))) %>% 
    mutate(across(everything(), ~ na_if(.x, "")))
}

v1v3_vsearch_tax <- import_qiime_taxonomy("~/R/nanopore_vs_illumina_2025/input/vsearch_zotuclassification_v1v3.qza")

v3v4_vsearch_tax <- import_qiime_taxonomy("~/R/nanopore_vs_illumina_2025/input/vsearch_zotuclassification_v3v4.qza")

ssumi_vsearch_tax <- import_qiime_taxonomy("~/R/nanopore_vs_illumina_2025/input/ssumi_vsearchZOTUclassification.qza")

# Sample metadata

metadata <- readxl::read_excel("~/R/nanopore_vs_illumina_2025/input/Nanopore_vs_Illumina_metadata.xlsx")

# Theoretical mock

theoretical_mock <- read_csv2("~/R/nanopore_vs_illumina_2025/input/Theoreticalmock2.csv") %>% 
  filter(Person == "skin_mock") %>%
  mutate(OTU = row_number())

theoretical_mock <- read_csv2("~/R/nanopore_vs_illumina_2025/input/Theoreticalmock2.csv") %>% 
  filter(Person == "staph_mock")%>%
  mutate(OTU = row_number())

skin_otu_theoretical <- theoretical_mock %>%
  select(OTU, Abundance) %>%
  column_to_rownames("OTU") %>%
  as.matrix()

skin_tax_theoretical <- theoretical_mock %>%
  select(OTU, Domain, Phylum, Class, Order, Family, Genus, Species) %>%
  column_to_rownames("OTU") %>%
  as.matrix()

staph_otu_theoretical <- theoretical_staph_mock %>%
  select(OTU, Abundance) %>%
  column_to_rownames("OTU") %>%
  as.matrix()

staph_tax_theoretical <- theoretical_staph_mock %>%
  select(OTU, Domain, Phylum, Class, Order, Family, Genus, Species) %>%
  column_to_rownames("OTU") %>%
  as.matrix()

# UMI counts

umi_counts <- read_tsv("~/R/nanopore_vs_illumina_2025/input/umi_counts.tsv") %>% 
  filter(Sample_Name != "barcode01")

############# REMOVE MITOCHONDRIA, CHLOROPLAST, UNKNOWN PHYLA ##################

v1v3_mc_otus <- v1v3_vsearch_tax %>% 
  filter(Family == "Mitochondria" | Order == "Chloroplast" | is.na(Phylum) | Domain == "Unassigned") %>% 
  pull(OTU)

v3v4_mc_otus <- v3v4_vsearch_tax %>% 
  filter(Family == "Mitochondria" | Order == "Chloroplast" | is.na(Phylum) | Domain == "Unassigned") %>%
  pull(OTU)

ssumi_mc_otus <- ssumi_vsearch_tax %>% 
  filter(Family == "Mitochondria" | Order == "Chloroplast" | is.na(Phylum) | Domain == "Unassigned") %>%
  pull(OTU)

v1v3_otutab <- v1v3_otutab[!rownames(v1v3_otutab) %in% v1v3_mc_otus, ]
v3v4_otutab <- v3v4_otutab[!rownames(v3v4_otutab) %in% v3v4_mc_otus, ]
ssumi_otutab <- ssumi_otutab[!rownames(ssumi_otutab) %in% ssumi_mc_otus, ]


######################## DECONTAMINATE PCR CONTROLS ############################

library(decotnam)

# Add a row of control = TRUE
metadata <- metadata %>% 
  mutate(is.neg = Person == "pcr_control")

physeq_v1v3 <- phyloseq(otu_table(as.matrix(v1v3_otutab), taxa_are_rows = TRUE),
                        sample_data(metadata %>% column_to_rownames("Sample"))) 

physeq_v3v4 <- phyloseq(otu_table(as.matrix(v3v4_otutab), taxa_are_rows = TRUE), 
                        sample_data(metadata %>% column_to_rownames("Sample"))) 

physeq_ssumi <- phyloseq(otu_table(as.matrix(ssumi_otutab), taxa_are_rows = TRUE), 
                         sample_data(metadata %>% column_to_rownames("Sample"))) 

contamdf_prev_v1v3 <- isContaminant(physeq_v1v3, method="prevalence", neg="is.neg", threshold=0.3)
table(contamdf_prev_v1v3$contaminant) 

contamdf_prev_v3v4 <- isContaminant(physeq_v3v4, method="prevalence", neg="is.neg", threshold=0.4) 
table(contamdf_prev_v3v4$contaminant) 

contamdf_prev_ssumi <- isContaminant(physeq_ssumi, method="prevalence", neg="is.neg", threshold=0.1) 
table(contamdf_prev_ssumi$contaminant) 


contaminants_v1v3 <- contamdf_prev_v1v3 %>%                                 
  filter(contaminant == TRUE) %>%
  rownames_to_column('name') %>%
  pull(name)

contaminants_v3v4 <- contamdf_prev_v3v4 %>%                             
  filter(contaminant == TRUE) %>%
  rownames_to_column('name') %>%
  pull(name)

contaminants_ssumi <- contamdf_prev_ssumi %>%                        
  filter(contaminant == TRUE) %>%
  rownames_to_column('name') %>%
  pull(name)

v1v3_filt_otutab <- v1v3_otutab
v1v3_filt_otutab[rownames(v1v3_filt_otutab) %in% contaminants_v1v3, ] <- 0

v3v4_filt_otutab <- v3v4_otutab
v3v4_filt_otutab[rownames(v3v4_filt_otutab) %in% contaminants_v3v4, ] <- 0

ssumi_filt_otutab <- ssumi_otutab
ssumi_filt_otutab[rownames(ssumi_filt_otutab) %in% contaminants_ssumi, ] <- 0

################### LIBRARY SIZE WITH "DECONTAM" METHOD ########################

v3v4_otutab_size <- v3v4_otutab %>% 
  t() %>%
  as.data.frame() %>% 
  rownames_to_column("Sample") %>%
  mutate(Library_Size = rowSums(across(where(is.numeric)))) %>%
  left_join(metadata, by = "Sample")


ggplot(v3v4_otutab_size, aes(x = Sample, y = Library_Size)) +
  geom_col() +
  labs(x = "Samples", y = "Library size") +
  theme_bw()

v1v3_otutab_size <- v1v3_filt_otutab %>% 
  t() %>%
  as.data.frame() %>% 
  rownames_to_column("Sample") %>%
  mutate(Library_Size = rowSums(across(where(is.numeric)))) %>%
  left_join(metadata, by = "Sample")

ggplot(v1v3_otutab_size, aes(x = Sample, y = Library_Size)) +
  geom_col() +
  labs(x = "Samples", y = "Library size") +
  theme_bw()

ssumi_otutab_size <- ssumi_filt_otutab %>% 
  t() %>%
  as.data.frame() %>% 
  rownames_to_column("Sample") %>%
  mutate(Library_Size = rowSums(across(where(is.numeric)))) %>%
  left_join(metadata, by = "Sample")

sum(ssumi_otutab_size$Library_Size)

ggplot(ssumi_otutab_size, aes(x = Sample, y = Library_Size)) +
  geom_col() +
  labs(x = "Samples", y = "Library size") +
  theme_bw()

############################# AVERAGE UMI COUNTS ###############################

ggplot(umi_counts, aes(x = ubs_count)) +
  theme_bw() +
  geom_bar(stat = "count") +
  facet_wrap(. ~ Sample_Name, ncol= 3, scales = "free_y") +
  labs(x = "UMI coverage", y = "Count")

umi_counts %>%
  filter(Sample_Name != "NK") %>%
  summarise(
    mean = mean(ubs_count),
    sd   = sd(ubs_count))

############################### RAREFACTION ####################################

v3v4_otutab_size %>%
  select(Sample, Library_Size) %>%
  arrange(desc(Library_Size))    

v1v3_otutab_size %>%
  select(Sample, Library_Size) %>%
  arrange(desc(Library_Size))

ssumi_otutab_size %>%
  select(Sample, Library_Size) %>%
  arrange(desc(Library_Size))

Srare_v3v4_zotu <- rarefy(v3v4_filt_otutab %>% t(), 179999)
Srrare_v3v4_zotu <- rrarefy(v3v4_filt_otutab %>% t(), 179999)
plot(specnumber(v3v4_filt_otutab %>% t()), 
     Srare_v3v4_zotu, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(v3v4_filt_otutab %>% t(), step = 200, sample = 179999, col = "blue", cex = 0.8, label = TRUE)

Srare_v1v3_zotu <- rarefy(v1v3_filt_otutab %>% t(), 179999)
Srrare_v1v3_zotu <- rrarefy(v1v3_filt_otutab %>% t(), 179999)
plot(specnumber(v1v3_filt_otutab %>% t()), 
     Srare_v1v3_zotu, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(v1v3_filt_otutab %>% t(), step = 200, sample = 179999, col = "blue", cex = 0.8, label = TRUE)

Srare_umi_zotu <- rarefy(ssumi_filt_otutab %>% t(), 15704)
Srrare_umi_zotu <- rrarefy(ssumi_filt_otutab %>% t(), 15704)
plot(specnumber(ssumi_filt_otutab %>% t()), 
     Srare_umi_zotu, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(ssumi_filt_otutab %>% t(), step = 50, sample = 15704, col = "blue", cex = 0.8, label = TRUE)

########################### PHYLOSEQ FROM RAREFIED ##############################

rarefied_physeq_v3v4 <- phyloseq(
  otu_table(Srrare_v3v4_zotu, taxa_are_rows = F),
  tax_table(v3v4_vsearch_tax %>% column_to_rownames("OTU") %>% as.matrix()),
  sample_data(metadata %>% column_to_rownames("Sample")))

rarefied_physeq_v1v3 <- phyloseq(
  otu_table(Srrare_v1v3_zotu, taxa_are_rows = F),
  tax_table(v1v3_vsearch_tax %>% column_to_rownames("OTU") %>% as.matrix()),
  sample_data(metadata %>% column_to_rownames("Sample")))

rarefied_physeq_ssumi <- phyloseq(
  otu_table(Srrare_umi_zotu, taxa_are_rows = F),
  tax_table(ssumi_vsearch_tax %>% column_to_rownames("OTU") %>% as.matrix()),
  sample_data(metadata %>% column_to_rownames("Sample")))

rarefied_physeq_v3v4_pruned <- prune_taxa(taxa_sums(rarefied_physeq_v3v4) > 0, rarefied_physeq_v3v4)
rarefied_physeq_v1v3_pruned <- prune_taxa(taxa_sums(rarefied_physeq_v1v3) > 0, rarefied_physeq_v1v3)
rarefied_physeq_ssumi_pruned <- prune_taxa(taxa_sums(rarefied_physeq_ssumi) > 0, rarefied_physeq_ssumi)

############################# ADD THEORETICAL MOCK #############################

# Skin mock phyloseq object
ps_skin_theoretical <- phyloseq(
  otu_table(skin_otu_theoretical, taxa_are_rows = TRUE),
  tax_table(skin_tax_theoretical))

# Staphylococcal mock phyloseq object
ps_staph_theoretical <- phyloseq(
  otu_table(staph_otu_theoretical, taxa_are_rows = TRUE),
  tax_table(staph_tax_theoretical))

################################################################################
#                           ALPHA DIVERSITY MEASURES                           #
################################################################################

################## THEORETICAL MOCK RICHNESS ON ZOTU LEVEL #####################

skin_richness_theoretical <- estimate_richness(ps_skin_theoretical, measures = NULL) %>%
  rownames_to_column("Sample") %>%
  mutate(
    Old_sample_name = "Theoretical",
    Primers = "Theoretical",
    Seq_method = "Theoretical",
    Person = "skin_mock",
    sample_mock_control = "skin_mock",
    Cycles = NA,
    Dilution = NA,
    is.neg = FALSE)

staph_richness_theoretical <- estimate_richness(ps_staph_theoretical, measures = NULL) %>%
  rownames_to_column("Sample") %>%
  mutate(
    Old_sample_name = "Theoretical",
    Primers = "Theoretical",
    Seq_method = "Theoretical",
    Person = "staph_mock",
    sample_mock_control = "staph_mock",
    Cycles = NA,
    Dilution = NA,
    is.neg = FALSE)

########################### RICHNESS ON ZOTU LEVEL #############################

richness_v3v4 <- estimate_richness(rarefied_physeq_v3v4_pruned, split = TRUE, measures = NULL) %>% 
  rownames_to_column("Sample") %>% 
  left_join(metadata)

richness_v1v3 <- estimate_richness(rarefied_physeq_v1v3_pruned, split = TRUE, measures = NULL) %>% 
  rownames_to_column(("Sample")) %>% 
  left_join(metadata)

richness_ssumi <- estimate_richness(rarefied_physeq_ssumi_pruned, split = TRUE, measures = NULL) %>% 
  rownames_to_column("Sample") %>% 
  left_join(metadata)

############ JOIN THEORETICAL AND SAMPLES RICHNESS ON ZOTU LEVEL ###############

alpha_div <- richness_v3v4 %>% 
  bind_rows(richness_v1v3) %>% 
  bind_rows(richness_ssumi) %>% 
  bind_rows(skin_richness_theoretical) %>%   
  bind_rows(staph_richness_theoretical)

alpha_div_plot <- alpha_div %>% 
  pivot_longer(cols = c("Observed", "Chao1", "se.chao1", "ACE", "se.ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"), names_to = "Index_name", values_to = "Index_value") %>% 
  mutate(Primers = factor(Primers, levels = c("Theoretical", "V1V9", "V1V3", "V3V4")))

##################### MOCK STATISTICS ON ZOTU LEVEL ############################

skin_alpha <- alpha_div_plot %>%
  filter(Person == "skin_mock") %>%
  filter(Index_name %in% c("Observed", "Shannon", "InvSimpson"))

staph_alpha <- alpha_div_plot %>%
  filter(Person == "staph_mock") %>%
  filter(Index_name %in% c("Observed", "Shannon", "InvSimpson"))

library(rstatix)

skin_alpha %>%
  filter(Index_name == "Shannon") %>%
  group_by(Primers) %>%
  summarise(
    mean = mean(Index_value),
    sd = sd(Index_value),
    n = n())

staph_alpha %>%
  filter(Index_name == "Shannon") %>%
  group_by(Primers) %>%
  summarise(
    mean = mean(Index_value),
    sd = sd(Index_value),
    n = n())

skin_ttest <- skin_alpha %>%
  group_by(Index_name) %>%
  pairwise_t_test(
    Index_value ~ Primers,
    p.adjust.method = "bonferroni",
    paired = FALSE,
    detailed = TRUE)

staph_ttest <- staph_alpha %>%
  group_by(Index_name) %>%
  pairwise_t_test(
    Index_value ~ Primers,
    p.adjust.method = "bonferroni",
    paired = FALSE,
    detailed = TRUE)

##################### SAMPLES STATISTICS ON ZOTU LEVEL #########################

sample_comparisons <- list(c("V1V9", "V1V3"), c("V1V9", "V3V4"), c("V1V3", "V3V4"))

samples_plot_paired <- alpha_div_plot %>%
  filter(Index_name %in% c("Observed", "Shannon", "InvSimpson")) %>%
  filter(sample_mock_control == "sample") %>%
  group_by(Person, Primers, Index_name) %>%
  summarise(Index_value = mean(Index_value), .groups = "drop") %>%
  mutate(Primers = factor(Primers, levels = c("V1V9", "V1V3", "V3V4"))) %>%
  
  ggplot(aes(x = Primers, y = Index_value, color = Primers)) +
  geom_boxplot() +
  geom_point(size = 2) +
  geom_line(aes(group = Person), color = "grey", linetype = "dashed") +
  facet_wrap(~ Index_name, scales = "free_y", ncol = 3) +
  stat_compare_means(comparisons = sample_comparisons,
                     method = "wilcox.test",
                     paired = TRUE,
                     label = "p.signif") +
  labs(title = "Sample alpha diversity comparison on ZOTU level", x = "Primer pair", y = "Alpha diversity index") +
  theme_bw() +
  theme(legend.position = "none")

############################### AGGLOMERATE GENUS ##############################

rarefied_physeq_v3v4_pruned_glom <- tax_glom(rarefied_physeq_v3v4_pruned, taxrank = "Genus", NArm = FALSE)
rarefied_physeq_v1v3_pruned_glom <- tax_glom(rarefied_physeq_v1v3_pruned, taxrank = "Genus", NArm = FALSE)
rarefied_physeq_ssumi_pruned_glom <- tax_glom(rarefied_physeq_ssumi_pruned, taxrank = "Genus", NArm = FALSE)
ps_staph_theoretical_glom <- tax_glom(ps_staph_theoretical, taxrank = "Genus", NArm = FALSE)
ps_skin_theoretical_glom <- tax_glom(ps_skin_theoretical, taxrank = "Genus", NArm = FALSE)

################################# RICHNESS #####################################

richness_v3v4_glom <- estimate_richness(rarefied_physeq_v3v4_pruned_glom, split = TRUE, measures = NULL) %>% 
  rownames_to_column("Sample") %>% 
  left_join(metadata)

richness_v1v3_glom <- estimate_richness(rarefied_physeq_v1v3_pruned_glom, split = TRUE, measures = NULL) %>% 
  rownames_to_column(("Sample")) %>% 
  left_join(metadata)

richness_ssumi_glom <- estimate_richness(rarefied_physeq_ssumi_pruned_glom, split = TRUE, measures = NULL) %>% 
  rownames_to_column("Sample") %>% 
  left_join(metadata)

##################### JOIN SAMPLES RICHNESS ON ZOTU LEVEL  #####################

alpha_div_glom <- richness_v3v4_glom %>% 
  bind_rows(richness_v1v3_glom) %>% 
  bind_rows(richness_ssumi_glom)

alpha_div_plot_glom <- alpha_div_glom %>% 
  pivot_longer(cols = c("Observed", "Chao1", "se.chao1", "ACE", "se.ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"), names_to = "Index_name", values_to = "Index_value") %>% 
  mutate(Primers = factor(Primers, levels = c("Theoretical", "V1V9", "V1V3", "V3V4")))

##################### STATISTICAL PLOT FOR SAMPLES #############################

sample_comparisons <- list(c("V1V9", "V1V3"), c("V1V9", "V3V4"), c("V1V3", "V3V4"))

samples_plot_paired_glom <- alpha_div_plot_glom %>%
  filter(Index_name %in% c("Observed", "Shannon", "InvSimpson")) %>%
  filter(sample_mock_control == "sample") %>%
  group_by(Person, Primers, Index_name) %>%
  summarise(Index_value = mean(Index_value), .groups = "drop") %>%
  mutate(Primers = factor(Primers, levels = c("V1V9", "V1V3", "V3V4"))) %>%
  
  ggplot(aes(x = Primers, y = Index_value, color = Primers)) +
  geom_boxplot() +
  geom_point(size = 2) +
  geom_line(aes(group = Person), color = "grey", linetype = "dashed") +
  facet_wrap(~ Index_name, scales = "free_y", ncol = 3) +
  stat_compare_means(comparisons = sample_comparisons,
                     method = "wilcox.test",
                     paired = TRUE,
                     label = "p.signif") +
  labs(title = "Sample alpha diversity comparison on genus level", x = "Primer pair", y = "Alpha diversity index") +
  theme_bw() +
  theme(legend.position = "none")

################################################################################
#                           BETA DIVERSITY MEASURES                            #
################################################################################

OTU_Srrare_v3v4_zotu <- as(otu_table(rarefied_physeq_v3v4_pruned), "matrix")
OTU_Srrare_v1v3_zotu <- as(otu_table(rarefied_physeq_v1v3_pruned), "matrix")
OTU_Srrare_ssumi_zotu <- as(otu_table(rarefied_physeq_ssumi_pruned), "matrix")

# Filter voulanteer samples

OTU_Srrare_v3v4_zotu_selected <- OTU_Srrare_v3v4_zotu[
  rownames(OTU_Srrare_v3v4_zotu) %in% c("NPCP_A1a", "NPCP_B1a2", "NPCP_C1a", "NPCP_D1a", "NPCP_E1a2"), , drop = FALSE]

OTU_Srrare_v1v3_zotu_selected <- OTU_Srrare_v1v3_zotu[
  rownames(OTU_Srrare_v1v3_zotu) %in% c("NPCP_A1b", "NPCP_B1b2", "NPCP_C1b", "NPCP_D1b", "NPCP_E1b2"), , drop = FALSE]

OTU_Srrare_umi_zotu_selected <- OTU_Srrare_ssumi_zotu[
  rownames(OTU_Srrare_ssumi_zotu) %in% c("MP_03_A1", "MP_03_B1", "MP_03_C1", "MP_03_D1", "MP_03_E1"), , drop = FALSE]

###################### DISSIMILARITY MATRICES ON ZOTU LEVEL ####################

set.seed(21)
dissimilarity_v3v4 <- avgdist(OTU_Srrare_v3v4_zotu_selected, 
                              sample = 179999, 
                              distfun = vegdist, 
                              meanfun = mean,
                              transf = NULL, 
                              iterations = 100, 
                              dmethod = "bray",
                              diag = TRUE, 
                              upper = FALSE)

set.seed(22)
dissimilarity_v1v3 <- avgdist(OTU_Srrare_v1v3_zotu_selected, 
                              sample = 179999, 
                              distfun = vegdist, 
                              meanfun = mean,
                              transf = NULL, 
                              iterations = 100, 
                              dmethod = "bray",
                              diag = TRUE, 
                              upper = FALSE)

set.seed(24)
dissimilarity_ssumi <- avgdist(OTU_Srrare_umi_zotu_selected, 
                               sample = 15704, 
                               distfun = vegdist, 
                               meanfun = mean,
                               transf = NULL, 
                               iterations = 100, 
                               dmethod = "bray",
                               diag = TRUE, 
                               upper = FALSE)

############################# PCoA ON ZOTU LEVEL ###############################

# V3-V4
pcoa_v3v4 <- cmdscale(dissimilarity_v3v4, eig = TRUE, add = TRUE)
positions_v3v4 <- pcoa_v3v4$points
colnames(positions_v3v4) <- c("PCoA1", "PCoA2")
percent_explained_v3v4 <- 100 * pcoa_v3v4$eig / sum(pcoa_v3v4$eig)
pretty_percent_explained_v3v4 <- format(round(percent_explained_v3v4, digits =1), nsmall=1, trim=TRUE)
labels_v3v4 <- c(glue("PCoA Axis 1 ({pretty_percent_explained_v3v4[1]}%)"),
                 glue("PCoA Axis 2 ({pretty_percent_explained_v3v4[2]}%)"))
positions_v3v4 <- positions_v3v4 %>% 
  as_tibble(rownames = "Sample") %>%
  left_join(metadata)

ggplot(positions_v3v4, aes(x = PCoA1, y = PCoA2)) + 
  geom_point(size = 4, ) +
  geom_text(aes(label = Sample)) +
  labs(x=labels_v3v4[1], y=labels_v3v4[2], title = "PCoA, Bray-Curtis") +
  stat_ellipse() +
  theme_bw()

# V1-V3
pcoa_v1v3 <- cmdscale(dissimilarity_v1v3, eig = TRUE, add = TRUE)
positions_v1v3 <- pcoa_v1v3$points
colnames(positions_v1v3) <- c("PCoA1", "PCoA2")
percent_explained_v1v3 <- 100 * pcoa_v1v3$eig / sum(pcoa_v1v3$eig)
pretty_percent_explained_v1v3 <- format(round(percent_explained_v1v3, digits =1), nsmall=1, trim=TRUE)
labels_v1v3 <- c(glue("PCoA Axis 1 ({pretty_percent_explained_v1v3[1]}%)"),
                 glue("PCoA Axis 2 ({pretty_percent_explained_v1v3[2]}%)"))
positions_v1v3 <- positions_v1v3 %>% 
  as_tibble(rownames = "Sample") %>%
  left_join(metadata)

ggplot(positions_v1v3, aes(x = PCoA1, y = PCoA2)) + 
  geom_point(size = 3, ) +
  geom_text(aes(label = Sample)) +
  labs(x=labels_v1v3[1], y=labels_v1v3[2], title = "PCoA, Bray-Curtis") +
  stat_ellipse() +
  theme_bw()

# ssUMI
pcoa_umi <- cmdscale(dissimilarity_ssumi, eig = TRUE, add = TRUE)
positions_umi <- pcoa_umi$points
colnames(positions_umi) <- c("PCoA1", "PCoA2")
percent_explained_umi <- 100 * pcoa_umi$eig / sum(pcoa_umi$eig)
pretty_percent_explained_umi <- format(round(percent_explained_umi, digits =1), nsmall=1, trim=TRUE)
labels_umi <- c(glue("PCoA Axis 1 ({pretty_percent_explained_umi[1]}%)"),
                glue("PCoA Axis 2 ({pretty_percent_explained_umi[2]}%)"))
positions_umi <- positions_umi %>% 
  as_tibble(rownames = "Sample") %>%
  left_join(metadata)

ggplot(positions_umi, aes(x = PCoA1, y = PCoA2)) + 
  geom_point(size = 3, ) +
  geom_text(aes(label = Sample)) +
  labs(x=labels_umi[1], y=labels_umi[2], title = "PCoA, Bray-Curtis") +
  stat_ellipse() +
  theme_bw()

########################### BETADISPER ON ZOTU LEVEL ###########################

# V3-V4
metadata_v3v4_filt <- metadata %>% 
  filter(Sample %in% c("NPCP_G1a", "NPCP_L1a2", "NPCP_M1a", "NPCP_N1a", "NPCP_U1a2"))

bd_v3v4 <- betadisper(dissimilarity_v3v4, metadata_v3v4_filt$sample_mock_control)

bd_dist_v3v4 <- bd_v3v4$distances %>%
  enframe(name = "Sample", value = "Dispersion")


# V1-V3
metadata_v1v3_filt <- metadata %>% 
  filter(Sample %in% c("NPCP_G1b", "NPCP_L1b2", "NPCP_M1b", "NPCP_N1b", "NPCP_U1b2"))

bd_v1v3 <- betadisper(dissimilarity_v1v3, metadata_v1v3_filt$sample_mock_control)

bd_dist_v1v3 <- bd_v1v3$distances %>%
  enframe(name = "Sample", value = "Dispersion")

# ssUMI
metadata_ssumi_filt <- metadata %>% 
  filter(Sample %in% c("MP_03_A1", "MP_03_B1", "MP_03_C1", "MP_03_D1", "MP_03_E1"))

bd_ssumi <- betadisper(dissimilarity_ssumi, metadata_ssumi_filt$sample_mock_control)

bd_dist_ssumi <- bd_ssumi$distances %>%
  enframe(name = "Sample", value = "Dispersion")

# Join
all_bd_distances <- bd_dist_v3v4 %>% 
  rbind(bd_dist_v1v3) %>% 
  rbind(bd_dist_ssumi) %>% 
  left_join(metadata)

# Plot
samples_beta_dispersion <- ggplot(all_bd_distances, 
                                  aes(x = factor(Primers, levels = c("V1V9", "V3V4", "V1V3")), 
                                      y = Dispersion)) +
  geom_line(aes(group = Person), color = "grey", linetype = "dashed") +
  geom_point(size = 4, aes(shape = Person)) +
  scale_shape_manual(values = c("A" = 2,
                                "B" = 1,
                                "C" = 4,
                                "D" = 3,
                                "E" = 5)) +
  theme_bw() +
  labs(y = "Distance to centroid on ZOTU level",
       x = "",
       title = "") +
  scale_x_discrete(labels = c("V1V9" = "ssUMI\nV1-V9\nONT",
                              "V3V4" = "V3-V4\nIllumina",
                              "V1V3" = "V1-V3\nIllumina")) +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 10, margin = margin(t = 10)),
        axis.title.y = element_text(size = 10, margin = margin(r = 10)),
        axis.text  = element_text(size = 8),
        panel.border = element_rect(linewidth = 0.5))

# Summarize
all_bd_distances_summary <- all_bd_distances %>% 
  group_by(Primers) %>%
  summarise(
    mean_disp = mean(Dispersion),
    sd_disp = sd(Dispersion),
    n = n())

write.table(all_bd_distances_summary, "~/R/nanopore_vs_illumina_2025/Paper_images_ideas/samples_beta_diversity_dispersion_on_zotu_lvl.csv")

########################## PROCRUSTES ON ZOTU LEVEL ############################

set.seed(126)
procord_v3v4_v1v3 <- procrustes(pcoa_v3v4, pcoa_v1v3, symmetric = TRUE)
set.seed(127)
procord_v1v3_ssumi <- procrustes(pcoa_v1v3, pcoa_umi, symmetric = TRUE)
set.seed(128)
procord_v3v4_ssumi  <- procrustes(pcoa_v3v4, pcoa_umi, symmetric = TRUE)

summary(procord_v3v4_v1v3)
summary(procord_v1v3_ssumi)
summary(procord_v3v4_ssumi)

protord_v3v4_v1v3 <- protest(pcoa_v3v4, pcoa_v1v3, symmetric = TRUE, permutations = 999)
protord_v1v3_ssumi <- protest(pcoa_v1v3, pcoa_umi, symmetric = TRUE, permutations = 999)
protord_v3v4_ssumi  <- protest(pcoa_v3v4, pcoa_umi, symmetric = TRUE, permutations = 999)

procrustes_ord_results <- data.frame(
  Comparison = c("V3V4 vs. V1V3", "V1V3 vs. ssUMI", "V3V4 vs. ssUMI"),
  m2_statistic = c(protord_v3v4_v1v3$ss,
                   protord_v1v3_ssumi$ss,
                   protord_v3v4_ssumi$ss),
  r_congruence = sqrt(1 - c(protord_v3v4_v1v3$ss,
                            protord_v1v3_ssumi$ss,
                            protord_v3v4_ssumi$ss)),
  p_value = c(protord_v3v4_v1v3$signif,
              protord_v1v3_ssumi$signif,
              protord_v3v4_ssumi$signif))

write.table(procrustes_ord_results, "~/R/nanopore_vs_illumina_2025/Paper_images_ideas/procrustes_ord_results_on_zotu_lvl.csv")

plot(procord_v3v4_v1v3, kind = 1, type = "text")
plot(procord_v1v3_ssumi, kind = 1, type = "text")
plot(procord_v3v4_ssumi, kind = 1, type = "text")

plot(procord_v3v4_v1v3, kind = 2)
plot(procord_v1v3_ssumi, kind = 2)
plot(procord_v3v4_ssumi, kind = 2)

######################### MANTEL TEST ON ZOTU LEVEL ############################

set.seed(129)
mantel_v3v4_v1v3 <- mantel(dissimilarity_v3v4, dissimilarity_v1v3,
                           method = "spearman", permutations = 999)
set.seed(130)
mantel_v1v3_ssumi <- mantel(dissimilarity_v1v3, dissimilarity_ssumi,
                            method = "spearman", permutations = 999)
set.seed(131)
mantel_v3v4_ssumi <- mantel(dissimilarity_v3v4, dissimilarity_ssumi,
                            method = "spearman", permutations = 999)

# Summary table
mantel_results <- data.frame(
  Comparison = c("V3V4 vs V1V3", "V1V3 vs ssUMI", "V3V4 vs ssUMI"),
  Mantel_r = c(mantel_v3v4_v1v3$statistic,
               mantel_v1v3_ssumi$statistic,
               mantel_v3v4_ssumi$statistic),
  p_value = c(mantel_v3v4_v1v3$signif,
              mantel_v1v3_ssumi$signif,
              mantel_v3v4_ssumi$signif))

write.table(mantel_results, "~/R/nanopore_vs_illumina_2025/Paper_images_ideas/mantel_test_results_on_zotu_lvl.csv")

################################# AGGLOMERATE TAXA ##############################

rarefied_physeq_v3v4_genus <- tax_glom(rarefied_physeq_v3v4_pruned,  taxrank = "Genus", NArm = FALSE)
rarefied_physeq_v1v3_genus <- tax_glom(rarefied_physeq_v1v3_pruned,  taxrank = "Genus", NArm = FALSE)
rarefied_physeq_ssumi_genus <- tax_glom(rarefied_physeq_ssumi_pruned, taxrank = "Genus", NArm = FALSE)

genus_Srrare_v3v4_zotu <- as(otu_table(rarefied_physeq_v3v4_genus), "matrix")
genus_Srrare_v1v3_zotu <- as(otu_table(rarefied_physeq_v1v3_genus), "matrix")
genus_Srrare_ssumi_zotu <- as(otu_table(rarefied_physeq_ssumi_genus), "matrix")

genus_Srrare_v3v4_zotu_selected <- genus_Srrare_v3v4_zotu[
  rownames(genus_Srrare_v3v4_zotu) %in% c("NPCP_A1a", "NPCP_B1a2", "NPCP_C1a", "NPCP_D1a", "NPCP_E1a2"), , drop = FALSE] 

genus_Srrare_v1v3_zotu_selected <- genus_Srrare_v1v3_zotu[
  rownames(genus_Srrare_v1v3_zotu) %in% c("NPCP_A1b", "NPCP_B1b2", "NPCP_C1b", "NPCP_D1b", "NPCP_E1b2"), , drop = FALSE]

genus_Srrare_umi_zotu_selected <- genus_Srrare_ssumi_zotu[
  rownames(genus_Srrare_ssumi_zotu) %in% c("MP_03_A1", "MP_03_B1", "MP_03_C1", "MP_03_D1", "MP_03_E1"), , drop = FALSE]

#################### DISSIMILARITY MATRICES ON GENUS LEVEL #####################

set.seed(21)
genus_dissimilarity_v3v4 <- avgdist(genus_Srrare_v3v4_zotu_selected, 
                                    sample = 179999, 
                                    distfun = vegdist, 
                                    meanfun = mean,
                                    transf = NULL, 
                                    iterations = 100, 
                                    dmethod = "bray",
                                    diag = TRUE, 
                                    upper = FALSE)

set.seed(22)
genus_dissimilarity_v1v3 <- avgdist(genus_Srrare_v1v3_zotu_selected, 
                                    sample = 179999, 
                                    distfun = vegdist, 
                                    meanfun = mean,
                                    transf = NULL, 
                                    iterations = 100, 
                                    dmethod = "bray",
                                    diag = TRUE, 
                                    upper = FALSE)

set.seed(24)
genus_dissimilarity_ssumi <- avgdist(genus_Srrare_umi_zotu_selected, 
                                     sample = 15704, 
                                     distfun = vegdist, 
                                     meanfun = mean,
                                     transf = NULL, 
                                     iterations = 100, 
                                     dmethod = "bray",
                                     diag = TRUE, 
                                     upper = FALSE)

############################# PCoA ON GENUS LEVEL ##############################

# V3-V4
genus_pcoa_v3v4 <- cmdscale(genus_dissimilarity_v3v4, eig = TRUE, add = TRUE)
genus_positions_v3v4 <- genus_pcoa_v3v4$points
colnames(genus_positions_v3v4) <- c("PCoA1", "PCoA2")
genus_percent_explained_v3v4 <- 100 * genus_pcoa_v3v4$eig / sum(genus_pcoa_v3v4$eig)
genus_pretty_percent_explained_v3v4 <- format(round(genus_percent_explained_v3v4, digits =1), nsmall=1, trim=TRUE)
genus_labels_v3v4 <- c(glue("PCoA Axis 1 ({genus_pretty_percent_explained_v3v4[1]}%)"),
                       glue("PCoA Axis 2 ({genus_pretty_percent_explained_v3v4[2]}%)"))
genus_positions_v3v4 <- genus_positions_v3v4 %>% 
  as_tibble(rownames = "Sample") %>%
  left_join(metadata)

ggplot(genus_positions_v3v4, aes(x = PCoA1, y = PCoA2)) + 
  geom_point(size = 4, ) +
  geom_text(aes(label = Sample)) +
  labs(x=labels_v3v4[1], y=labels_v3v4[2], title = "PCoA, Bray-Curtis") +
  stat_ellipse() +
  theme_bw()

# V1-V3
genus_pcoa_v1v3 <- cmdscale(genus_dissimilarity_v1v3, eig = TRUE, add = TRUE)
genus_positions_v1v3 <- genus_pcoa_v1v3$points
colnames(genus_positions_v1v3) <- c("PCoA1", "PCoA2")
genus_percent_explained_v1v3 <- 100 * genus_pcoa_v1v3$eig / sum(genus_pcoa_v1v3$eig)
genus_pretty_percent_explained_v1v3 <- format(round(genus_percent_explained_v1v3, digits =1), nsmall=1, trim=TRUE)
genus_labels_v1v3 <- c(glue("PCoA Axis 1 ({genus_pretty_percent_explained_v1v3[1]}%)"),
                       glue("PCoA Axis 2 ({genus_pretty_percent_explained_v1v3[2]}%)"))
genus_positions_v1v3 <- genus_positions_v1v3 %>% 
  as_tibble(rownames = "Sample") %>%
  left_join(metadata)

ggplot(genus_positions_v1v3, aes(x = PCoA1, y = PCoA2)) + 
  geom_point(size = 4, ) +
  geom_text(aes(label = Sample)) +
  labs(x=labels_v1v3[1], y=labels_v1v3[2], title = "PCoA, Bray-Curtis") +
  stat_ellipse() +
  theme_bw()

# ssUMI
genus_pcoa_umi <- cmdscale(genus_dissimilarity_ssumi, eig = TRUE, add = TRUE)
genus_positions_umi <- genus_pcoa_umi$points
colnames(genus_positions_umi) <- c("PCoA1", "PCoA2")
genus_percent_explained_umi <- 100 * genus_pcoa_umi$eig / sum(genus_pcoa_umi$eig)
genus_pretty_percent_explained_umi <- format(round(genus_percent_explained_umi, digits =1), nsmall=1, trim=TRUE)
genus_labels_umi <- c(glue("PCoA Axis 1 ({genus_pretty_percent_explained_umi[1]}%)"),
                      glue("PCoA Axis 2 ({genus_pretty_percent_explained_umi[2]}%)"))
genus_positions_umi <- genus_positions_umi %>% 
  as_tibble(rownames = "Sample") %>%
  left_join(metadata)

ggplot(genus_positions_umi, aes(x = PCoA1, y = PCoA2)) + 
  geom_point(size = 4, ) +
  geom_text(aes(label = Sample)) +
  labs(x=labels_umi[1], y=labels_umi[2], title = "PCoA, Bray-Curtis") +
  stat_ellipse() +
  theme_bw()

############################ BETADISPER ON GENUS LEVEL #########################

# V3-V4
metadata_v3v4_filt <- metadata %>% 
  filter(Sample %in% c("NPCP_G1a", "NPCP_L1a2", "NPCP_M1a", "NPCP_N1a", "NPCP_U1a2"))

genus_bd_v3v4 <- betadisper(genus_dissimilarity_v3v4, metadata_v3v4_filt$sample_mock_control)

genus_bd_dist_v3v4 <- genus_bd_v3v4$distances %>%
  enframe(name = "Sample", value = "Dispersion")

# V1-V3
metadata_v1v3_filt <- metadata %>% 
  filter(Sample %in% c("NPCP_G1b", "NPCP_L1b2", "NPCP_M1b", "NPCP_N1b", "NPCP_U1b2"))

genus_bd_v1v3 <- betadisper(genus_dissimilarity_v1v3, metadata_v1v3_filt$sample_mock_control)

genus_bd_dist_v1v3 <- genus_bd_v1v3$distances %>%
  enframe(name = "Sample", value = "Dispersion")

# ssUMI
metadata_ssumi_filt <- metadata %>% 
  filter(Sample %in% c("MP_03_A1", "MP_03_B1", "MP_03_C1", "MP_03_D1", "MP_03_E1"))

genus_bd_ssumi <- betadisper(genus_dissimilarity_ssumi, metadata_ssumi_filt$sample_mock_control)

genus_bd_dist_ssumi <- genus_bd_ssumi$distances %>%
  enframe(name = "Sample", value = "Dispersion")

genus_all_bd_distances <- genus_bd_dist_v3v4 %>% 
  rbind(genus_bd_dist_v1v3) %>% 
  rbind(genus_bd_dist_ssumi) %>% 
  left_join(metadata)

genus_samples_beta_dispersion <- ggplot(genus_all_bd_distances, 
                                        aes(x = factor(Primers, levels = c("V1V9", "V3V4", "V1V3")), 
                                            y = Dispersion)) +
  geom_line(aes(group = Person), color = "grey", linetype = "dashed") +
  geom_point(size = 4, aes(shape = Person)) +
  scale_shape_manual(values = c("A" = 2,
                                "B" = 1,
                                "C" = 4,
                                "D" = 3,
                                "E" = 5)) +
  theme_bw() +
  labs(y = "Distance to centroid on genus level",
       x = "",
       title = "") +
  scale_x_discrete(labels = c("V1V9" = "ssUMI\nV1-V9\nONT",
                              "V3V4" = "V3-V4\nIllumina",
                              "V1V3" = "V1-V3\nIllumina")) +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 10, margin = margin(t = 10)),
        axis.title.y = element_text(size = 10, margin = margin(r = 10)),
        axis.text  = element_text(size = 8),
        panel.border = element_rect(linewidth = 0.5))

combined_plot <- cowplot::plot_grid(samples_beta_dispersion, genus_samples_beta_dispersion, ncol = 2, labels = c("A", "B"))

ggsave("~/R/nanopore_vs_illumina_2025/Paper_images_ideas/beta_dispersion_combined.pdf", combined_plot, width = 170, height = 100, units = "mm", dpi = 600)
ggsave("~/R/nanopore_vs_illumina_2025/Paper_images_ideas/beta_dispersion_combined.svg", combined_plot, width = 170, height = 100, units = "mm", dpi = 600)

genus_all_bd_distances_summary <- genus_all_bd_distances %>% 
  group_by(Primers) %>%
  summarise(
    mean_disp = mean(Dispersion),
    sd_disp = sd(Dispersion),
    n = n())

write.table(genus_all_bd_distances_summary, "~/R/nanopore_vs_illumina_2025/Paper_images_ideas/samples_beta_diversity_dispersion_on_genus_lvl.csv")

######################### PROCRUSTES ON GENUS LEVEL ############################

genus_pcoa_v3v4
genus_pcoa_v1v3
genus_pcoa_umi

# Procrustes — with ordiantions
set.seed(126)
genus_procord_v3v4_v1v3 <- procrustes(genus_pcoa_v3v4, pcoa_v1v3, symmetric = TRUE)
set.seed(127)
genus_procord_v1v3_ssumi <- procrustes(genus_pcoa_v1v3, pcoa_umi, symmetric = TRUE)
set.seed(128)
genus_procord_v3v4_ssumi  <- procrustes(genus_pcoa_v3v4, pcoa_umi, symmetric = TRUE)

summary(genus_procord_v3v4_v1v3)
summary(genus_procord_v1v3_ssumi)
summary(genus_procord_v3v4_ssumi)

# Protest for significance testing
genus_protord_v3v4_v1v3 <- protest(genus_pcoa_v3v4, pcoa_v1v3, symmetric = TRUE, permutations = 999)
genus_protord_v1v3_ssumi <- protest(genus_pcoa_v1v3, pcoa_umi, symmetric = TRUE, permutations = 999)
genus_protord_v3v4_ssumi  <- protest(genus_pcoa_v3v4, pcoa_umi, symmetric = TRUE, permutations = 999)


# Summary
genus_procrustes_ord_results <- data.frame(
  Comparison = c("V3V4 vs. V1V3", "V1V3 vs. ssUMI", "V3V4 vs. ssUMI"),
  m2_statistic = c(genus_protord_v3v4_v1v3$ss,
                   genus_protord_v1v3_ssumi$ss,
                   genus_protord_v3v4_ssumi$ss),
  r_congruence = sqrt(1 - c(genus_protord_v3v4_v1v3$ss,
                            genus_protord_v1v3_ssumi$ss,
                            genus_protord_v3v4_ssumi$ss)),
  p_value = c(genus_protord_v3v4_v1v3$signif,
              genus_protord_v1v3_ssumi$signif,
              genus_protord_v3v4_ssumi$signif))

write.table(genus_procrustes_ord_results, "~/R/nanopore_vs_illumina_2025/Paper_images_ideas/procrustes_ord_results_on_genus_lvl.csv")

plot(genus_procord_v3v4_v1v3, kind = 1, type = "text")
plot(genus_procord_v1v3_ssumi, kind = 1, type = "text")
plot(genus_procord_v3v4_ssumi, kind = 1, type = "text")

plot(genus_procord_v3v4_v1v3, kind = 2)
plot(genus_procord_v1v3_ssumi, kind = 2)
plot(genus_procord_v3v4_ssumi, kind = 2)

######################## MANTEL TEST ON GENUS LEVEL ############################

set.seed(129)
genus_mantel_v3v4_v1v3 <- mantel(genus_dissimilarity_v3v4, dissimilarity_v1v3,
                                 method = "spearman", permutations = 999)
set.seed(130)
genus_mantel_v1v3_ssumi <- mantel(genus_dissimilarity_v1v3, dissimilarity_ssumi,
                                  method = "spearman", permutations = 999)
set.seed(131)
genus_mantel_v3v4_ssumi <- mantel(genus_dissimilarity_v3v4, dissimilarity_ssumi,
                                  method = "spearman", permutations = 999)

# Summary table
genus_mantel_results <- data.frame(
  Comparison = c("V3V4 vs V1V3", "V1V3 vs ssUMI", "V3V4 vs ssUMI"),
  Mantel_r= c(genus_mantel_v3v4_v1v3$statistic,
              genus_mantel_v1v3_ssumi$statistic,
              genus_mantel_v3v4_ssumi$statistic),
  p_value = c(genus_mantel_v3v4_v1v3$signif,
              genus_mantel_v1v3_ssumi$signif,
              genus_mantel_v3v4_ssumi$signif))

print(genus_mantel_results)

write.table(genus_mantel_results, "~/R/nanopore_vs_illumina_2025/Paper_images_ideas/mantel_test_results_on_genus_lvl.csv")

################################################################################
#                   SKIN MOCK COMMUNITY COMPOSITION ON SPECIES LEVEL           #
################################################################################

species_in_skin_mock <- theoretical_mock %>% 
  select(Species)

genus_in_skin_mock <- theoretical_mock %>% 
  select(Genus)

taxa_names(rarefied_physeq_v3v4)  <- paste0("v3v4_",  taxa_names(rarefied_physeq_v3v4))
taxa_names(rarefied_physeq_v1v3)  <- paste0("v1v3_",  taxa_names(rarefied_physeq_v1v3))
taxa_names(rarefied_physeq_ssumi) <- paste0("ssUMI_", taxa_names(rarefied_physeq_ssumi))

all_primers_phyloseq <- merge_phyloseq(rarefied_physeq_v3v4,
                                       rarefied_physeq_v1v3,
                                       rarefied_physeq_ssumi)

all_primers_phyloseq_skin_mock <- all_primers_phyloseq %>% 
  subset_samples(Person == "skin_mock")

melt_all_primers_phyloseq_skin_mock <- psmelt(all_primers_phyloseq_skin_mock)

theoretical_plus_melt <- melt_all_primers_phyloseq_skin_mock %>% 
  bind_rows(theoretical_mock)

####################### SKIN MOCK COMMUNITY TAXA PLOT ##########################

# Select species that appear in the mock community, agglomerate other species (NA, sp., metagenome, uncultured, ...)
theoretical_plus_melt <- theoretical_plus_melt %>% 
  mutate(Species2 = case_when(
    Species %in% species_in_skin_mock$Species ~ Species,
    is.na(Species) ~ "Other",
    TRUE ~ "Other")) %>% 
  mutate(Genus2 = case_when(
    Genus %in% genus_in_skin_mock$Genus ~ Genus,
    is.na(Genus) ~ "Other",
    TRUE ~ "Other"))

how_many_other <- theoretical_plus_melt %>%
  group_by(Primers, Species2) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Primers) %>%
  mutate(Percentage = round(Abundance / sum(Abundance) * 100, 2)) %>%
  arrange(Primers, desc(Percentage))  

# Remove other genera for plotting
theoretical_plus_melt <- theoretical_plus_melt %>% 
  filter(Genus2 != "Other")

taxon_colors <- c(
  # Cutibacterium
  "Cutibacterium acnes" = "#054A7A",
  "Cutibacterium avidum" = "#098BD9",
  "Cutibacterium (other)" = "#56B4E9",
  # Corynebacterium
  "Corynebacterium striatum" = "#6BAED6",
  "Corynebacterium afermentas" = "#C6DBEF",
  "Corynebacterium (other)" = "#9ECAE1",
  # Micrococcus
  "Micrococcus (other)" = "#9E9AC8",
  "Micrococcus luteus" = "#DADAEB",
  # Kocuria
  "Kocuria marina" = "#CC79A7",
  "Kocuria (other)" = "#EFB6D6",
  # Dermacoccus
  "Dermacoccus nishinomiyaensis" = "#807DBA",
  "Dermacoccus (other)" = "#6A51A3",
  # Acinetobacter
  "Acinetobacter baumannii" = "#148F77",
  "Acinetobacter bereziniae" = "#43BA8F",
  "Acinetobacter (other)" = "#A3E4D7",
  # Pseudomonas
  "Pseudomonas putida" = "#74C476",
  "Pseudomonas aeruginosa" = "#C7E9C0",
  "Pseudomonas (other)" = "#238B45",
  # Roseomonas
  "Roseomonas mucosa" = "#BDB0AE",
  "Roseomonas (other)" = "#DBD6D5",
  # Staphylococcus
  "Staphylococcus aureus" = "#FF7F00",
  "Staphylococcus capitis" = "#F09163",
  "Staphylococcus epidermidis" = "#FE9929",
  "Staphylococcus haemolyticus" = "#FDAE6B",
  "Staphylococcus hominis" = "#FEC44F",
  "Staphylococcus succinus" = "#FFD5AF",
  "Staphylococcus warneri" = "#FEEDA0",
  "Staphylococcus (other)" = "#9D654C",
  # Deinococcus
  "Deinococcus radiodurans" = "#c9645d",
  "Deinococcus (other)" = "#964c47")

theoretical_plus_melt <- theoretical_plus_melt %>%
  mutate(Species2 = gsub("_", " ", Species2))

theoretical_plus_melt <- theoretical_plus_melt %>%
  mutate(Species2 = if_else(
    Species2 == "Other",
    paste0(Genus2, " (other)"),
    Species2))

theoretical_plus_melt$Primers <- factor(theoretical_plus_melt$Primers, 
                                        levels = c("Theoretical", "V1V9", "V3V4", "V1V3"))
theoretical_plus_melt$Seq_method2 <- factor(theoretical_plus_melt$Seq_method, 
                                            levels = c("Theoretical", "Nanopore", "Illumina"))

skin_mock_composition <- 
  ggplot(theoretical_plus_melt, aes(x = Primers, y = Abundance, color = Species2, fill = Species2)) +
  geom_col(position = "fill") +
  scale_x_discrete(labels = c("V1V9" = "ssUMI\nV1-V9\nONT",
                              "V3V4" = "V3-V4\nIllumina",
                              "V1V3" = "V1-V3\nIllumina")) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = taxon_colors) +
  scale_color_manual(values = taxon_colors, guide = "none") +
  labs(y = "Relative abundance (%)", x = "", fill = "Species in the skin mock community") +
  theme_bw() +
  theme(
    legend.position = "bottom", 
    legend.key.size = unit(0.25, "cm"), 
    legend.key.spacing.x = unit(0, "pt"),
    legend.key.spacing.y = unit(0, "pt"),
    legend.title.position = "top",
    legend.box = "vertical",
    axis.title = element_text(size = 9),
    axis.text  = element_text(size = 8),
    legend.text = element_text(size = 6.3),
    legend.title = element_text(size = 9),
    panel.border = element_rect(linewidth = 0.5)) +
  guides(fill = guide_legend(byrow = FALSE, ncol = 2))

################ SKIN MOCK COMMUNITY L1 NORM ON SPECIES LEVEL ##################

l1_plot_data_replicates <- theoretical_plus_melt %>%
  mutate(Species2 = gsub("_", " ", Species2)) %>%
  filter(Primers != "Theoretical") %>%  # remove theoretical
  group_by(Primers, Species2, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Primers, Sample) %>%
  mutate(rel_abund = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  select(-Abundance) %>% 
  mutate(Replicate = case_when(
    Sample %in% c("Skin_MOCK_1", "NPCP_SKIN1b", "NPCP_SKIN1a") ~ "Rep1",
    Sample %in% c("Skin_MOCK_2", "NPCP_SKIN2b", "NPCP_SKIN2a") ~ "Rep2",
    Sample %in% c("Skin_MOCK_3", "NPCP_SKIN3b", "NPCP_SKIN3a") ~ "Rep3"))

l1_plot_data_wide <- l1_plot_data_replicates %>%
  select(-Sample) %>%
  pivot_wider(names_from = Primers,
              values_from = rel_abund,
              values_fill = 0)

# Average L1
l1_values <- l1_plot_data_wide %>%
  group_by(Replicate) %>%
  summarise(
    L1_V1V9_V1V3 = sum(abs(V1V9 - V1V3)),
    L1_V1V9_V3V4 = sum(abs(V1V9 - V3V4)),
    L1_V1V3_V3V4 = sum(abs(V1V3 - V3V4))) %>%
  summarise(across(starts_with("L1"), mean))

# L1 for each species
l1_plot_data_avg <- l1_plot_data_wide %>%
  group_by(Species2) %>%
  summarise(
    V1V9 = mean(V1V9),
    V1V3 = mean(V1V3),
    V3V4 = mean(V3V4))

################# SKIN MOCK COMMUNITY L1 NORM ON GENUS LEVEL ###################

# data wrangle
l1_genus_replicates <- theoretical_plus_melt %>%
  mutate(Genus2 = gsub("_", " ", Genus2)) %>%
  filter(Primers != "Theoretical") %>%
  group_by(Primers, Genus2, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Primers, Sample) %>%
  mutate(rel_abund = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  select(-Abundance) %>%
  mutate(Replicate = case_when(
    Sample %in% c("Skin_MOCK_1", "NPCP_SKIN1b", "NPCP_SKIN1a") ~ "Rep1",
    Sample %in% c("Skin_MOCK_2", "NPCP_SKIN2b", "NPCP_SKIN2a") ~ "Rep2",
    Sample %in% c("Skin_MOCK_3", "NPCP_SKIN3b", "NPCP_SKIN3a") ~ "Rep3"))

l1_genus_wide <- l1_genus_replicates %>%
  select(-Sample) %>%
  pivot_wider(names_from = Primers,
              values_from = rel_abund,
              values_fill = 0)

# Average L1
l1_genus_values <- l1_genus_wide %>%
  group_by(Replicate) %>%
  summarise(
    L1_V1V9_V1V3 = sum(abs(V1V9 - V1V3)),
    L1_V1V9_V3V4 = sum(abs(V1V9 - V3V4)),
    L1_V1V3_V3V4 = sum(abs(V1V3 - V3V4))
  ) %>%
  summarise(across(starts_with("L1"), mean))

# L1 for each genus
l1_genus_data_avg <- l1_genus_wide %>%
  group_by(Genus2) %>%
  summarise(
    V1V9 = mean(V1V9),
    V1V3 = mean(V1V3),
    V3V4 = mean(V3V4))

genus_colors <- c(
  "Cutibacterium" = "#054A7A",
  "Corynebacterium" = "#9ECAE1",
  "Micrococcus" = "#9E9AC8",
  "Kocuria" = "#CC79A7",
  "Dermacoccus" = "#6A51A3",
  "Acinetobacter" = "#148F77",
  "Pseudomonas" = "#C7E9C0",
  "Roseomonas" = "#DBD6D5",
  "Staphylococcus" = "#FE9929")

# Plot first comparison
p1_genus <- ggplot(l1_genus_data_avg, aes(x = V1V9, y = V1V3, label = Genus2)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  geom_point(shape = 21, size = 4, color = "black", aes(fill = Genus2)) +
  scale_fill_manual(values = genus_colors, name = "Genus") +  
  annotate("text", x = 0, y = 0.5 * 0.95, label = paste0("L1 = ", round(l1_genus_values$L1_V1V9_V1V3, 2)), hjust = 0, size = 3) +
  labs(x = "ssUMI V1-V9 ONT", y = "V1-V3 Illumina", title = "") +
  xlim(0, 0.5) + 
  ylim(0, 0.5) +
  theme_bw()+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text  = element_text(size = 6),
        panel.border = element_rect(linewidth = 0.5)) +
  coord_fixed(ratio = 1)

# Plot second comparison
p2_genus <- ggplot(l1_genus_data_avg, aes(x = V1V9, y = V3V4, label = Genus2)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  geom_point(shape = 21, size = 4, color = "black", aes(fill = Genus2)) +
  scale_fill_manual(values = genus_colors, name = "Genus") + 
  annotate("text", x = 0, y = 0.5 * 0.95, label = paste0("L1 = ", round(l1_genus_values$L1_V1V9_V3V4, 2)), hjust = 0, size = 3) +
  labs(x = "ssUMI V1-V9 ONT", y = "V3-V4 Illumina", title = "") +
  xlim(0, 0.5) + 
  ylim(0, 0.5) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text  = element_text(size = 6),
        panel.border = element_rect(linewidth = 0.5)) +
  coord_fixed(ratio = 1)

# Plot third comparison (this one keeps its legend, which we'll extract)
p3_genus <- ggplot(l1_genus_data_avg, aes(x = V1V3, y = V3V4, label = Genus2)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  geom_point(shape = 21, size = 4, color = "black", aes(fill = Genus2)) +
  scale_fill_manual(values = genus_colors, name = "Genus") + 
  annotate("text", x = 0, y = 0.5 * 0.95, label = paste0("L1 = ", round(l1_genus_values$L1_V1V3_V3V4, 2)), hjust = 0, size = 3) +
  labs(x = "V1-V3 Illumina", y = "V3-V4 Illumina", title = "") +
  theme_bw()+
  xlim(0, 0.5) + 
  ylim(0, 0.5) +
  theme(legend.position = "right",
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text  = element_text(size = 6),
        panel.border = element_rect(linewidth = 0.5)) +
  coord_fixed(ratio = 1)

# Join
top_row <- cowplot::plot_grid(p1_genus, p2_genus, ncol = 2, labels = c("A", "B"))
bottom_row <- cowplot::plot_grid(p3_genus, ncol = 1, labels = c("C"))

combined_l1_genus <- cowplot::plot_grid(top_row, bottom_row, ncol = 1)

ggsave("~/R/nanopore_vs_illumina_2025/Paper_images_ideas/l1_genus_combined.pdf", combined_l1_genus, width = 170, height = 170, units = "mm", dpi = 600)
ggsave("~/R/nanopore_vs_illumina_2025/Paper_images_ideas/l1_genus_combined.svg", combined_l1_genus, width = 170, height = 170, units = "mm", dpi = 600)

####### AVERAGE WEIGHTED CLASSIFICATION CONSENSUS IN SKIN MOCK COMMUNITY #######

theoretical_plus_melt_weighted <- theoretical_plus_melt %>% 
  mutate(Consensus = as.numeric(Consensus)) %>% 
  filter(Abundance > 0, !is.na(Consensus)) %>% 
  group_by(Primers, Species2) %>% 
  summarise(
    n_rows = n_distinct(OTU),
    weighted_consensus = weighted.mean(Consensus, w = Abundance),
    unweighted_consensus = mean(Consensus),
    .groups = "drop")

# Plot aspect ratio
n_species_skin <- length(unique(theoretical_plus_melt_weighted$Species2[theoretical_plus_melt_weighted$Primers != "Theoretical"]))
n_primers <- 3  # V1V3, V3V4, V1V9

# Plot

library(ggh4x)

all_skin_species <- theoretical_plus_melt %>%
  distinct(Species2) %>%
  pull(Species2)

heatmap1 <- theoretical_plus_melt_weighted %>% 
  filter(Primers != "Theoretical", Species2 != "Other (other)") %>%
  complete(Primers, Species2 = all_skin_species) %>%  # expand to full species list
  mutate(Species2 = factor(Species2, levels = rev(sort(unique(Species2))))) %>% 
  
  ggplot(aes(x = factor(Primers, levels = c("V1V9", "V3V4", "V1V3")), y = Species2, fill = weighted_consensus)) +
  geom_tile(width = 0.95, height = 0.95, color = "white", linewidth = 0.3) +
  scale_fill_gradient(low = "#F2F2F6", high = "#26267C",
                      name = "Weighted consensus",
                      limits = c(0.5, 1),
                      na.value = "transparent") +
  scale_x_discrete(labels = c("V1V9" = "ssUMI\nV1-V9 ONT",
                              "V3V4" = "V3-V4\nIllumina",
                              "V1V3" = "V1-V3\nIllumina")) +
  theme_bw() +
  labs(x = "", y = "") +
  theme(
    legend.position = "none",
    panel.border = element_rect(linewidth = 0.5)) +
  force_panelsizes(rows = unit(length(all_skin_species) * 4, "mm"),
                   cols = unit(n_primers * 16, "mm"))

######### SKIN MOCK COMMUNITY ALPHA DIVERISTY MEASURES ON SPECIES LEVEL ########

mock_grouped_complete_wider <- theoretical_plus_melt %>% 
  ungroup() %>% 
  filter(!is.na(Sample)) %>% 
  select(Species2, Abundance, Sample) %>% 
  pivot_wider(names_from = Sample, values_from = Abundance, values_fn = list(Abundance = sum), values_fill = list(Abundance = 0)) %>% 
  mutate(across(where(is.numeric), ~replace_na(.x, 0)))

mock_grouped_complete_wider <- mock_grouped_complete_wider %>% 
  column_to_rownames(var = "Species2")

podatki <- theoretical_plus_melt %>%
  ungroup() %>%                   
  select(Sample, Primers) %>%
  distinct() %>% 
  filter(!is.na(Sample))

alpha_mock <- mock_grouped_complete_wider %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  pivot_longer(cols=-Sample) %>% 
  group_by(Sample) %>%
  summarize(richness = specnumber(value),
            shannon= diversity(value, index = "shannon"),
            simpson = diversity(value, index = "simpson"),
            inversesimpson = 1/simpson,
            n = sum(value))

alpha_mock <- alpha_mock %>% 
  left_join(podatki, by = join_by(Sample))

shannon_ttest <- alpha_mock %>%
  pairwise_t_test(
    shannon ~ Primers,
    p.adjust.method = "bonferroni",
    paired = FALSE,
    detailed = TRUE)

################################################################################
#              STAPHYLOCOCCAL MOCK COMMUNITY COMPOSITION ON SPECIES LEVEL      #
################################################################################

species_in_staph_mock <- theoretical_staph_mock %>% 
  select(Species)

genus_in_staph_mock <- theoretical_staph_mock %>% 
  select(Genus)

taxa_names(rarefied_physeq_v3v4)  <- paste0("v3v4_",  taxa_names(rarefied_physeq_v3v4))
taxa_names(rarefied_physeq_v1v3)  <- paste0("v1v3_",  taxa_names(rarefied_physeq_v1v3))
taxa_names(rarefied_physeq_ssumi) <- paste0("ssUMI_", taxa_names(rarefied_physeq_ssumi))

all_primers_phyloseq <- merge_phyloseq(rarefied_physeq_v3v4,
                                       rarefied_physeq_v1v3,
                                       rarefied_physeq_ssumi)

all_primers_phyloseq_staph_mock <- all_primers_phyloseq %>% 
  subset_samples(Person == "staph_mock")

melt_all_primers_phyloseq_staph_mock <- psmelt(all_primers_phyloseq_staph_mock)

staph_theoretical_plus_melt <- melt_all_primers_phyloseq_staph_mock %>% 
  bind_rows(theoretical_staph_mock)

################### STAPHYLOCOCCAL MOCK COMMUNITY TAXA PLOT ####################

# Select species that appear in the mock community, agglomerate other species (NA, sp., metagenome, uncultured, ...)
taph_theoretical_plus_melt <- staph_theoretical_plus_melt %>% 
  mutate(Species2 = case_when(
    Species %in% species_in_staph_mock$Species ~ Species,      
    is.na(Species) ~ "Other",                             
    TRUE ~ "Other"              
  )) %>% 
  mutate(Genus2 = case_when(
    Genus %in% genus_in_staph_mock$Genus ~ Genus,
    is.na(Genus) ~ "Other",
    TRUE ~ "Other"))

how_many_other_staph <- staph_theoretical_plus_melt %>%
  group_by(Primers, Genus2, Species2) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Primers) %>%
  mutate(Percentage = round(Abundance / sum(Abundance) * 100, 2)) %>%
  arrange(Primers, desc(Percentage))

# Remove other genera for plotting
staph_theoretical_plus_melt <- staph_theoretical_plus_melt %>% 
  filter(Genus2 != "Other")

staph_taxon_colors <- c(
  "Staphylococcus aureus"        = "#FF7F00",
  "Staphylococcus capitis"       = "#F09163",
  "Staphylococcus caprae"        = "#D94701",
  "Staphylococcus epidermidis"   = "#FE9929",
  "Staphylococcus haemolyticus"  = "#FDAE6B",
  "Staphylococcus hominis"       = "#FEC44F",
  "Staphylococcus lugdunensis"   = "#FFF0C0",
  "Staphylococcus pasteuri"      = "#FFE873",
  "Staphylococcus petrasii"      = "#F2ce00",
  "Staphylococcus succinus"      = "#FFD5AF",
  "Staphylococcus warneri"       = "#FEEDA0",
  "Staphylococcus (other)"       = "#9D654C")

staph_theoretical_plus_melt <- staph_theoretical_plus_melt %>%
  mutate(Species2 = gsub("_", " ", Species2))

staph_theoretical_plus_melt <- staph_theoretical_plus_melt %>%
  mutate(Species2 = if_else(
    Species2 == "Other",
    paste0(Genus2, " (other)"),
    Species2))

staph_mock_composition <-
  ggplot(staph_theoretical_plus_melt, aes(x = Primers, y = Abundance, color = Species2, fill = Species2)) +
  geom_col(position = "fill") +
  scale_x_discrete(labels = c("V1V9" = "ssUMI\nV1-V9\nONT",
                              "V3V4" = "V3-V4\nIllumina",
                              "V1V3" = "V1-V3\nIllumina")) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = staph_taxon_colors) +
  scale_color_manual(values = staph_taxon_colors, guide = "none") +
  labs(y = "Relative abundance (%)", x = "", fill = "Species in the staphylococcal mock community") +
  theme_bw() +
  theme(
    legend.position = "bottom", 
    legend.key.size = unit(0.25, "cm"), 
    legend.key.spacing.x = unit(0, "pt"),
    legend.key.spacing.y = unit(0, "pt"),
    legend.title.position = "top",
    legend.box = "vertical",
    axis.title = element_text(size = 9),
    axis.text  = element_text(size = 8),
    legend.text = element_text(size = 6.3),
    legend.title = element_text(size = 9),
    panel.border = element_rect(linewidth = 0.5)) +
  guides(fill = guide_legend(byrow = FALSE, ncol = 2))

# Join skin and staphylococcal mocks

combined_mock <- cowplot::plot_grid(skin_mock_composition, 
                                    staph_mock_composition,
                                    ncol = 2, 
                                    labels = c("A", "B"),
                                    align = "h")  # aligns plot panels horizontally

ggsave("~/R/nanopore_vs_illumina_2025/Paper_images_ideas/mock_composition_side_by_side.pdf", combined_mock, width = 170, height = 150, units = "mm", dpi = 600)
ggsave("~/R/nanopore_vs_illumina_2025/Paper_images_ideas/mock_composition_side_by_side.svg", combined_mock, width = 170, height = 150, units = "mm", dpi = 600)

########## STAPHYLOCOCCAL MOCK COMMUNITY L1 NORM ON SPECIES LEVEL ##############

staph_l1_plot_data_replicates <- staph_theoretical_plus_melt %>%
  mutate(Species2 = gsub("_", " ", Species2)) %>%
  filter(Primers != "Theoretical") %>%  # remove theoretical
  group_by(Primers, Species2, Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Primers, Sample) %>%
  mutate(rel_abund = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  select(-Abundance) %>% 
  mutate(Replicate = case_when(
    Sample %in% c("Staph_MOCK_1", "NPCP_STAPH1b", "NPCP_STAPH1a") ~ "Rep1",
    Sample %in% c("Staph_MOCK_2", "NPCP_STAPH2b", "NPCP_STAPH2a") ~ "Rep2",
    Sample %in% c("Staph_MOCK_3", "NPCP_STAPH3b", "NPCP_STAPH3a") ~ "Rep3"))

staph_l1_plot_data_wide <- staph_l1_plot_data_replicates %>%
  select(-Sample) %>%
  pivot_wider(names_from = Primers,
              values_from = rel_abund,
              values_fill = 0)

# Average L1
staph_l1_values <- staph_l1_plot_data_wide %>%
  group_by(Replicate) %>%
  summarise(
    L1_V1V9_V1V3 = sum(abs(V1V9 - V1V3)),
    L1_V1V9_V3V4 = sum(abs(V1V9 - V3V4)),
    L1_V1V3_V3V4 = sum(abs(V1V3 - V3V4))) %>%
  summarise(across(starts_with("L1"), mean)) 

# L1 for each species
staph_l1_plot_data_avg <- staph_l1_plot_data_wide %>%
  group_by(Species2) %>%
  summarise(
    V1V9 = mean(V1V9),
    V1V3 = mean(V1V3),
    V3V4 = mean(V3V4))

## AVERAGE WEIGHTED CLASSIFICATION CONSENSUS IN STAPHYLOCOCCAL MOCK COMMUNITY ##

staph_theoretical_plus_melt_weighted <- staph_theoretical_plus_melt %>% 
  mutate(Consensus = as.numeric(Consensus)) %>% 
  filter(Abundance > 0, !is.na(Consensus)) %>% 
  group_by(Primers, Species2) %>% 
  summarise(
    n_rows = n_distinct(OTU),
    weighted_consensus = weighted.mean(Consensus, w = Abundance),
    unweighted_consensus = mean(Consensus),
    .groups = "drop")

# Plot aspect ratio
n_species_staph <- length(unique(staph_theoretical_plus_melt_weighted$Species2[staph_theoretical_plus_melt_weighted$Primers != "Theoretical"]))
n_primers <- 3  # V1V3, V3V4, V1V9

# Plot

all_staph_species <- staph_theoretical_plus_melt %>%
  distinct(Species2) %>%
  pull(Species2)

heatmap2 <- staph_theoretical_plus_melt_weighted %>% 
  filter(Primers != "Theoretical", Species2 != "Other (other)", !is.na(Primers)) %>%
  mutate(Primers = droplevels(factor(Primers))) %>%
  complete(Primers, Species2 = all_staph_species) %>% # expand to full species list
  mutate(Species2 = factor(Species2, levels = rev(sort(unique(Species2))))) %>% 
  
  ggplot(aes(x = factor(Primers, levels = c("V1V9", "V3V4", "V1V3")), y = Species2, fill = weighted_consensus)) +
  geom_tile(width = 0.95, height = 0.95, color = "white", linewidth = 0.3) +
  scale_fill_gradient(low = "#F2F2F6", high = "#26267C",
                      name = "Weighted consensus",
                      limits = c(0.5, 1),
                      na.value = "transparent") +
  scale_x_discrete(labels = c("V1V9" = "ssUMI\nV1-V9 ONT",
                              "V3V4" = "V3-V4\nIllumina",
                              "V1V3" = "V1-V3\nIllumina")) +
  theme_bw() +
  labs(x = "", y = "") +
  theme(
    legend.position = "bottom",
    legend.title.position = "top",
    panel.border = element_rect(linewidth = 0.5)) +
  force_panelsizes(rows = unit(length(all_staph_species) * 4, "mm"),
                   cols = unit(n_primers * 16, "mm"))

# Join
n_rows_top <- length(all_skin_species)
n_rows_bottom <- length(all_staph_species)

combined_heatmaps <- cowplot::plot_grid(heatmap1, heatmap2, ncol = 1, align = "v", labels = c("A", "B"), rel_heights = c(n_rows_top, n_rows_bottom + 5))

ggsave("~/R/nanopore_vs_illumina_2025/Paper_images_ideas/combined_heatmaps.pdf", combined_heatmaps, width = 170, height = 250, units = "mm", dpi = 600)
ggsave("~/R/nanopore_vs_illumina_2025/Paper_images_ideas/combined_heatmaps.svg", combined_heatmaps, width = 170, height = 250, units = "mm", dpi = 600)

#### STAPHYLOCOCCAL MOCK COMMUNITY ALPHA DIVERISTY MEASURES ON SPECIES LEVEL ###

staph_mock_grouped_complete_wider <- staph_theoretical_plus_melt %>% 
  ungroup() %>% 
  filter(!is.na(Sample)) %>% 
  select(Species2, Abundance, Sample) %>% 
  pivot_wider(names_from = Sample, values_from = Abundance, values_fn = list(Abundance = sum), values_fill = list(Abundance = 0)) %>% 
  mutate(across(where(is.numeric), ~replace_na(.x, 0)))

staph_mock_grouped_complete_wider <- staph_mock_grouped_complete_wider %>% 
  column_to_rownames(var = "Species2")

staph_podatki <- staph_theoretical_plus_melt %>%
  ungroup() %>%                    # remove all previous grouping
  select(Sample, Primers) %>%
  distinct() %>% 
  filter(!is.na(Sample))

staph_alpha_mock <- staph_mock_grouped_complete_wider %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  pivot_longer(cols=-Sample) %>% 
  group_by(Sample) %>%
  summarize(richness = specnumber(value),
            shannon= diversity(value, index = "shannon"),
            simpson = diversity(value, index = "simpson"),
            inversesimpson = 1/simpson,
            n = sum(value))

staph_alpha_mock <- staph_alpha_mock %>% 
  left_join(staph_podatki, by = join_by(Sample))

staph_shannon_ttest <- staph_alpha_mokc %>%
  pairwise_t_test(
    shannon ~ Primers,
    p.adjust.method = "bonferroni",
    paired = FALSE,
    detailed = TRUE)



