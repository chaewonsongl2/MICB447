#### R SCRIPT - HISEAS DATASET - ALPHA DIVERSITY PLOTS

# Load CRAN packages
library(tidyverse)
library(vegan)
library(ape)

# Load Bioconductor packages
library(phyloseq)

# Load additional ggplot packages
library(ggplot2)
library(ggthemes)

# Define the set of random numbers
set.seed(800)

# Define name of metadata file
metadata_file <- "hiseas_metadata.txt"
# Load the metadata file
dat <- read_tsv(metadata_file)
# Convert time_sampling_days from a continuous variable to new categorical variable (resupply_event) with 3 groups
dat <- mutate(dat, resupply_event = cut(time_sampling_days,
                                                   breaks = c(0, 28, 70, Inf),
                                                   labels = c("Event 1 (day 15)","Event 2 (day 43)","Event 3 (day 107)")))
# Save new metadata file with new category
write_tsv(dat, "hiseas_metadata_resupply_event.txt")

# Import data
resupply_biom_file <- import_biom("resupply-table-with-taxonomy.biom")
plastic_wood_biom_file <- import_biom("plastic-wood-table-with-taxonomy.biom")
shower_biom_file <- import_biom("shower-table-with-taxonomy.biom")
metadata <- import_qiime_sample_data("hiseas_metadata_resupply_event.txt")
tree <- read_tree_greengenes("tree.nwk")

# Change time_sampling_days to a categorical variable
metadata$time_sampling_days <- as.factor(metadata$time_sampling_days)

# Convert from multichotomous to dichotomous tree
tree <- multi2di(tree) 

# Gather taxonomy, metadata, and phylogenetic tree information into phyloseq objects
resupply_physeq <- merge_phyloseq(resupply_biom_file, metadata, tree)
plastic_wood_physeq <- merge_phyloseq(plastic_wood_biom_file, metadata, tree)
shower_physeq <- merge_phyloseq(shower_biom_file, metadata, tree)

# Generate two shower physeq objects filtered on crew members 32 and 34 and remove exact shower timing samples
shower_32_physeq <- subset_samples(shower_physeq, crew_id == "32")
shower_32_physeq <- subset_samples(shower_32_physeq, showers_timing != "exact")
shower_34_physeq <- subset_samples(shower_physeq, crew_id == "34")
shower_34_physeq <- subset_samples(shower_34_physeq, showers_timing != "exact")

# Rarefy to the rarefaction depth
resupply_physeq_rar <- rarefy_even_depth(resupply_physeq, sample.size = 20706)
plastic_wood_physeq_rar <- rarefy_even_depth(plastic_wood_physeq, sample.size = 30000)
shower_32_physeq_rar <- rarefy_even_depth(shower_32_physeq, sample.size = 40000)
shower_34_physeq_rar <- rarefy_even_depth(shower_34_physeq, sample.size = 40000)

#### PLOT ALPHA DIVERSITY METRICS FOR EACH SAMPLE IN EACH RARIFIED FILTERED DATASET

## For resupply events

# Compute Shannon diversity of samples from each sampling day 
resupply_shannon <- plot_richness(resupply_physeq_rar, x="time_sampling_days", measures = c("Shannon"))

# Convert the Shannon metrics computed for resupply events into a data frame
resupply_shannon_table <- resupply_shannon$data

# Generate boxplots of Shannon diversity with fill color representing whether the sample was obtained before or after a resupply event
resupply_shannon_plot <- ggplot(resupply_shannon_table, aes(x = time_sampling_days, y = value, fill = supply_timing)) +
  geom_boxplot() +
  labs(title = "Shannon diversity on abiotic surfaces before and after 3 resupply events",
       x     = "Sampling day",
       y     = "Shannon Diversity") + 
  guides(fill = guide_legend("Supply timing")) +
  facet_wrap(~ resupply_event, scales = "free_x") +
  theme_few() +
  theme(axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13), axis.title = element_text(size = 16), legend.text = element_text(size = 15), legend.title = element_text(size = 16), strip.text = element_text(size = 15))

ggsave("resupply_shannon_plot.png", resupply_shannon_plot, width = 10, height = 6)

## For surface material 

# Compute Shannon diversity of samples from each surface material 
material_shannon <- plot_richness(plastic_wood_physeq_rar, x="orig_env_material", measures = c("Shannon")) 

# Convert the Shannon metrics computed for each surface material into a data frame
material_shannon_table <- material_shannon$data

# Plot Shannon diversity of samples from each surface material with fill color representing whether the sample was obtained from a wood or plastic surface
material_shannon_plot <- ggplot(material_shannon_table, aes(x = orig_env_material, y = value, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Shannon diversity on different abiotic surface materials",
       x     = "Abiotic surface material",
       y     = "Shannon Diversity") + 
  guides(fill = FALSE) +
  theme_few() +
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 13), axis.title = element_text(size = 17))

ggsave("material_shannon_plot.png", material_shannon_plot, width = 6, height = 6)

## For shower timing - crew member 32

# Compute Shannon diversity of skin samples from crew member 32 obtained before or after showering
shower_32_shannon <- plot_richness(shower_32_physeq_rar, x="showers_timing", measures = c("Shannon")) 

# Convert the Shannon metrics computed for each skin sample into a data frame
shower_32_shannon_table <- shower_32_shannon$data

# Plot Shannon diversity of skin samples from crew member 32 obtained before or after showering with fill colour representing shower timing
shower_32_shannon_plot <- ggplot(shower_32_shannon_table, aes(x = showers_timing, y = value, fill = showers_timing)) +
  geom_boxplot() +
  scale_x_discrete(limits=c("before", "after")) +
  labs(title = "Shannon diversity on skin of crew member 32 before and after showering",
       x     = "Shower timing",
       y     = "Shannon Diversity") + 
  guides(fill = FALSE) +
  theme_few() +
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 13), axis.title = element_text(size = 17))

ggsave("shower_crew32_shannon_plot.png", shower_32_shannon_plot, width = 6, height = 6)

## For shower timing - crew member 34

# Compute Shannon diversity of skin samples from crew member 34 obtained before or after showering
shower_34_shannon <- plot_richness(shower_34_physeq_rar, x="showers_timing", measures = c("Shannon")) 

# Convert the Shannon metrics computed for each skin sample into a data frame
shower_34_shannon_table <- shower_34_shannon$data

# Plot Shannon diversity of skin samples from crew member 32 obtained before or after showering with fill colour representing shower timing
shower_34_shannon_plot <- ggplot(shower_34_shannon_table, aes(x = showers_timing, y = value, fill = showers_timing)) +
  geom_boxplot() +
  scale_x_discrete(limits=c("before", "after")) +
  labs(title = "Shannon diversity on skin of crew member 34 before and after showering",
       x     = "Shower timing",
       y     = "Shannon Diversity") + 
  guides(fill = FALSE) +
  theme_few() +
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 13), axis.title = element_text(size = 17))

ggsave("shower_crew34_shannon_plot.png", shower_34_shannon_plot, width = 6, height = 6)




