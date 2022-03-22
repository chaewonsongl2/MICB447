#### R SCRIPT - HISEAS DATASET - BETA DIVERSITY PLOTS

# Load CRAN packages
library(tidyverse)
library(vegan)
library(ape)

# Load Bioconductor packages
library(phyloseq)
library(DESeq2)

# Load additional ggplot packages
library(ggplot2)
library(ggthemes)

# Define the set of random numbers
set.seed(800)

# Import data
resupply_biom_file <- import_biom("resupply-table-with-taxonomy.biom")
plastic_wood_biom_file <- import_biom("plastic-wood-table-with-taxonomy.biom")
shower_biom_file <- import_biom("shower-table-with-taxonomy.biom")
metadata <- import_qiime_sample_data("hiseas_metadata2.txt")
tree <- read_tree_greengenes("tree.nwk")

# Convert the multichotomous tree to a dichotomous tree
tree <- multi2di(tree)

# Gather taxonomy, metadata, and phylogenetic tree information into phyloseq objects
resupply_physeq <- merge_phyloseq(resupply_biom_file, metadata, tree)
plastic_wood_physeq <- merge_phyloseq(plastic_wood_biom_file, metadata, tree)
shower_physeq <- merge_phyloseq(shower_biom_file, metadata, tree)

# Assign names to numeric taxonomic ranks 
colnames(tax_table(resupply_physeq)) <- c("Kingdom", "Phylum", "Class","Order", "Family", "Genus", "Species")
colnames(tax_table(plastic_wood_physeq)) <- c("Kingdom", "Phylum", "Class","Order", "Family", "Genus", "Species")
colnames(tax_table(shower_physeq)) <- c("Kingdom", "Phylum", "Class","Order", "Family", "Genus", "Species")

# Filter shower physeq object by crew members 32 and 34 
shower_32_physeq <- subset_samples(shower_physeq, crew_id == "32")
shower_32_physeq <- subset_samples(shower_32_physeq, showers_timing != "exact" )
shower_34_physeq <- subset_samples(shower_physeq, crew_id == "34")
shower_34_physeq <- subset_samples(shower_34_physeq, showers_timing != "exact" )

# Rarefy to the the rarefaction depth 
resupply_physeq_rar <- rarefy_even_depth(resupply_physeq, sample.size = 20706)
plastic_wood_physeq_rar <- rarefy_even_depth(plastic_wood_physeq, sample.size = 30000)
shower_32_physeq_rar <- rarefy_even_depth(shower_32_physeq, sample.size = 40000)
shower_34_physeq_rar <- rarefy_even_depth(shower_34_physeq, sample.size = 40000)

#### PERFORM BETA DIVERSITY ANALYSIS WITH ORDINATE

# For resupply events
resupply_wunifrac_ord <- ordinate(resupply_physeq_rar, method = "PCoA", distance = "wunifrac")
resupply_uunifrac_ord <- ordinate(resupply_physeq_rar, method = "PCoA", distance = "uunifrac")
resupply_bray_ord <- ordinate(resupply_physeq_rar, method = "PCoA", distance = "bray")
resupply_jaccard_ord <- ordinate(resupply_physeq_rar, method = "PCoA", distance = "jaccard")

# For material 
plastic_wood_wunifrac_ord <- ordinate(plastic_wood_physeq_rar, method = "PCoA", distance = "wunifrac")
plastic_wood_uunifrac_ord <- ordinate(plastic_wood_physeq_rar, method = "PCoA", distance = "uunifrac")
plastic_wood_bray_ord <- ordinate(plastic_wood_physeq_rar, method = "PCoA", distance = "bray")
plastic_wood_jaccard_ord <- ordinate(plastic_wood_physeq_rar, method = "PCoA", distance = "jaccard")

# For shower timing - crew member 32
shower_32_wunifrac_ord <- ordinate(shower_32_physeq_rar, method = "PCoA", distance = "wunifrac")
shower_32_uunifrac_ord <- ordinate(shower_32_physeq_rar, method = "PCoA", distance = "uunifrac")
shower_32_bray_ord <- ordinate(shower_32_physeq_rar, method = "PCoA", distance = "bray")
shower_32_jaccard_ord <- ordinate(shower_32_physeq_rar, method = "PCoA", distance = "jaccard")

# For shower timing - crew member 34
shower_34_wunifrac_ord <- ordinate(shower_34_physeq_rar, method = "PCoA", distance = "wunifrac")
shower_34_uunifrac_ord <- ordinate(shower_34_physeq_rar, method = "PCoA", distance = "uunifrac")
shower_34_bray_ord <- ordinate(shower_34_physeq_rar, method = "PCoA", distance = "bray")
shower_34_jaccard_ord <- ordinate(shower_34_physeq_rar, method = "PCoA", distance = "jaccard")

#### GENERATE PCOA PLOTS

# For resupply events
pcoa_resupply_wunifrac <- plot_ordination(resupply_physeq_rar,
                resupply_wunifrac_ord,
                color = "time_sampling_days", shape = "sample_origin") +
  labs(title = "Resupply Event PCoA (Weighted UniFrac)") +
  theme_few() +
  guides(color = guide_legend("Sampling Day"), shape = guide_legend("Sampling Location")) +
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11), axis.title = element_text(size = 12), legend.text = element_text(size = 13), legend.title = element_text(size = 13))

ggsave("pcoa_resupply_wunifrac.png", pcoa_resupply_wunifrac, width = 8, height = 6)

pcoa_resupply_uunifrac <- plot_ordination(resupply_physeq_rar,
                resupply_uunifrac_ord,
                color = "time_sampling_days", shape = "sample_origin") +
  labs(title = "Resupply Event PCoA (Unweighted UniFrac)") +
  scale_color_manual(breaks = c("14", "28", "42", "56", "98", "115"),
                     values=c("blue4", "blue1", "dodgerblue3", "skyblue3", "skyblue2", "lightskyblue1")) +
  theme_few() +
  guides(color = guide_legend("Sampling Day"), shape = guide_legend("Sampling Location")) +
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11), axis.title = element_text(size = 12), legend.text = element_text(size = 13), legend.title = element_text(size = 13))

ggsave("pcoa_resupply_uunifrac.png", pcoa_resupply_uunifrac, width = 8, height = 6)

pcoa_resupply_bray <- plot_ordination(resupply_physeq_rar,
                resupply_bray_ord,
                color = "time_sampling_days", shape = "sample_origin") +
  labs(title = "Resupply Event PCoA (Bray-Curtis)") +
  theme_few() +
  guides(color = guide_legend("Sampling Day"), shape = guide_legend("Sampling Location")) +
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11), axis.title = element_text(size = 12), legend.text = element_text(size = 13), legend.title = element_text(size = 13))

ggsave("pcoa_resupply_bray.png", pcoa_resupply_bray, width = 8, height = 6)

pcoa_resupply_jaccard <- plot_ordination(resupply_physeq_rar,
                resupply_jaccard_ord,
                color = "time_sampling_days", shape = "sample_origin") +
  labs(title = "Resupply Event PCoA (Jaccard)") +
  theme_few() +
  guides(color = guide_legend("Sampling Day"), shape = guide_legend("Sampling Location")) +
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11), axis.title = element_text(size = 12), legend.text = element_text(size = 13), legend.title = element_text(size = 13))

ggsave("pcoa_resupply_jaccard.png", pcoa_resupply_jaccard, width = 8, height = 6)

# For surface material
pcoa_plastic_wood_wunifrac <- plot_ordination(plastic_wood_physeq_rar,
                plastic_wood_wunifrac_ord,
                color = "orig_env_material") +
  labs(title = "Surface Material PCoA (weighted UniFrac)") +
  theme_few() +
  stat_ellipse(type = "norm", size = 0.5) +
  guides(color = guide_legend("Surface Material")) +
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11), axis.title = element_text(size = 12), legend.text = element_text(size = 13), legend.title = element_text(size = 13))

ggsave("pcoa_plastic_wood_wunifrac.png", pcoa_plastic_wood_wunifrac, width = 8, height = 6)

pcoa_plastic_wood_uunifrac <- plot_ordination(plastic_wood_physeq_rar,
                plastic_wood_uunifrac_ord,
                color = "orig_env_material") +
  labs(title = "Surface Material PCoA (unweighted UniFrac)") +
  theme_few() +
  stat_ellipse(type = "norm", size = 0.5) +
  guides(color = guide_legend("Surface Material")) +
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11), axis.title = element_text(size = 12), legend.text = element_text(size = 13), legend.title = element_text(size = 13))

ggsave("pcoa_plastic_wood_uunifrac.png", pcoa_plastic_wood_uunifrac, width = 8, height = 6)

pcoa_plastic_wood_bray <- plot_ordination(plastic_wood_physeq_rar,
                plastic_wood_bray_ord,
                color = "orig_env_material") +
  labs(title = "Surface Material PCoA (Bray-Curtis)") +
  theme_few() +
  stat_ellipse(type = "norm", size = 0.5) +
  guides(color = guide_legend("Surface Material")) +
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11), axis.title = element_text(size = 12), legend.text = element_text(size = 13), legend.title = element_text(size = 13))

ggsave("pcoa_plastic_wood_bray.png", pcoa_plastic_wood_bray, width = 8, height = 6)

pcoa_plastic_wood_jaccard <- plot_ordination(plastic_wood_physeq_rar,
                plastic_wood_jaccard_ord,
                color = "orig_env_material") +
  labs(title = "Surface Material PCoA (Jaccard)") +
  theme_few() +
  stat_ellipse(type = "norm", size = 0.5) +
  guides(color = guide_legend("Surface Material")) +
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11), axis.title = element_text(size = 12), legend.text = element_text(size = 13), legend.title = element_text(size = 13))

ggsave("pcoa_plastic_wood_jaccard.png", pcoa_plastic_wood_jaccard, width = 8, height = 6)

# For shower timing - crew member 32
pcoa_shower_32_wunifrac <- plot_ordination(shower_32_physeq_rar,
                shower_32_wunifrac_ord,
                color = "showers_timing") +
  labs(title = "Shower Timing PCoA (weighted UniFrac) Crew 32") +
  theme_few() +
  guides(color = guide_legend("Shower Timing")) +
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11), axis.title = element_text(size = 12), legend.text = element_text(size = 13), legend.title = element_text(size = 13))

ggsave("pcoa_shower_32_wunifrac.png", pcoa_shower_32_wunifrac, width = 8, height = 6)

pcoa_shower_32_uunifrac <- plot_ordination(shower_32_physeq_rar,
                shower_32_uunifrac_ord,
                color = "showers_timing") +
  labs(title = "Shower Timing PCoA (unweighted UniFrac) Crew 32") +
  theme_few() +
  guides(color = guide_legend("Shower Timing")) +
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11), axis.title = element_text(size = 12), legend.text = element_text(size = 13), legend.title = element_text(size = 13))

ggsave("pcoa_shower_32_uunifrac.png", pcoa_shower_32_uunifrac, width = 8, height = 6)

pcoa_shower_32_bray <- plot_ordination(shower_32_physeq_rar,
                shower_32_bray_ord,
                color = "showers_timing") +
  labs(title = "Shower Timing PCoA (Bray-Curtis) Crew 32") +
  theme_few() +
  guides(color = guide_legend("Shower Timing")) +
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11), axis.title = element_text(size = 12), legend.text = element_text(size = 13), legend.title = element_text(size = 13))

ggsave("pcoa_shower_32_bray.png", pcoa_shower_32_bray, width = 8, height = 6)

pcoa_shower_32_jaccard <- plot_ordination(shower_32_physeq_rar,
                shower_32_jaccard_ord,
                color = "showers_timing") +
  labs(title = "Shower Timing PCoA (Jaccard) Crew 32") +
  theme_few() +
  guides(color = guide_legend("Shower Timing")) +
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11), axis.title = element_text(size = 12), legend.text = element_text(size = 13), legend.title = element_text(size = 13))

ggsave("pcoa_shower_32_jaccard.png", pcoa_shower_32_jaccard, width = 8, height = 6)

# For crew member 32 - PCoA with unweighted unifrac, shapes = shower timing, colour = time_sampling_days
pcoa_shower_32_uunifrac_sampling_day <- plot_ordination(shower_32_physeq_rar,
                                           shower_32_uunifrac_ord,
                                           color = "time_sampling_days", shape = "showers_timing") +
  labs(title = "Shower Timing PCoA (unweighted UniFrac) Crew 32") +
  theme_few() +
  guides(color = guide_legend("Sampling Day"), shape = guide_legend("Shower Timing")) +
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11), axis.title = element_text(size = 12), legend.text = element_text(size = 13), legend.title = element_text(size = 13))

ggsave("pcoa_shower_32_uunifrac_sampling_day.png", pcoa_shower_32_uunifrac_sampling_day, width = 8, height = 6)

# For shower timing - crew member 34
pcoa_shower_34_wunifrac <- plot_ordination(shower_34_physeq_rar,
                                           shower_34_wunifrac_ord,
                                           color = "showers_timing") +
  labs(title = "Shower Timing PCoA (weighted UniFrac) Crew 34") +
  theme_few() +
  guides(color = guide_legend("Shower Timing")) +
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11), axis.title = element_text(size = 12), legend.text = element_text(size = 13), legend.title = element_text(size = 13))

ggsave("pcoa_shower_34_wunifrac.png", pcoa_shower_34_wunifrac, width = 8, height = 6)

pcoa_shower_34_uunifrac <- plot_ordination(shower_34_physeq_rar,
                                           shower_34_uunifrac_ord,
                                           color = "showers_timing") +
  labs(title = "Shower Timing PCoA (unweighted UniFrac) Crew 34") +
  theme_few() +
  guides(color = guide_legend("Shower Timing")) +
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11), axis.title = element_text(size = 12), legend.text = element_text(size = 13), legend.title = element_text(size = 13))

ggsave("pcoa_shower_34_uunifrac.png", pcoa_shower_34_uunifrac, width = 8, height = 6)

pcoa_shower_34_bray <- plot_ordination(shower_34_physeq_rar,
                                       shower_34_bray_ord,
                                       color = "showers_timing") +
  labs(title = "Shower Timing PCoA (Bray-Curtis) Crew 34") +
  theme_few() +
  guides(color = guide_legend("Shower Timing")) +
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11), axis.title = element_text(size = 12), legend.text = element_text(size = 13), legend.title = element_text(size = 13))

ggsave("pcoa_shower_34_bray.png", pcoa_shower_34_bray, width = 8, height = 6)

pcoa_shower_34_jaccard <- plot_ordination(shower_34_physeq_rar,
                                          shower_34_jaccard_ord,
                                          color = "showers_timing") +
  labs(title = "Shower Timing PCoA (Jaccard) Crew 34") +
  theme_few() +
  guides(color = guide_legend("Shower Timing")) +
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11), axis.title = element_text(size = 12), legend.text = element_text(size = 13), legend.title = element_text(size = 13))

ggsave("pcoa_shower_34_jaccard.png", pcoa_shower_34_jaccard, width = 8, height = 6)

# For crew member 34 - PCoA with unweighted unifrac, shapes = shower timing, colour = time_sampling_days
pcoa_shower_34_uunifrac_sampling_day <- plot_ordination(shower_34_physeq_rar,
                                                        shower_34_uunifrac_ord,
                                                        color = "time_sampling_days", shape = "showers_timing") +
  labs(title = "Shower Timing PCoA (unweighted UniFrac) Crew 34") +
  theme_few() +
  guides(color = guide_legend("Sampling Day"), shape = guide_legend("Shower Timing")) +
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11), axis.title = element_text(size = 12), legend.text = element_text(size = 13), legend.title = element_text(size = 13))

ggsave("pcoa_shower_34_uunifrac_sampling_day.png", pcoa_shower_34_uunifrac_sampling_day, width = 8, height = 6)
