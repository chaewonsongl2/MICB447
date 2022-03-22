#### QIIME2 SHELL SCRIPT - HISEAS DATASET

# Create a directory in the data folder and navigate to it
mkdir /data/hiseas_analysis
cd /data/hiseas_analysis

# Import data using the manifest format while working directory is "/data/hiseas_analysis"
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path /mnt/datasets/project_2/hiseas/hiseas_manifest.txt \
  --output-path ./demux_seqs.qza
  
# Summarize the demultiplexed sequences to generate a visualization file
# Input the visualization file into QIIME2 View and determine an optimal truncation length from the quality plot
qiime demux summarize \
  --i-data ./demux_seqs.qza \
  --o-visualization ./demux_seqs.qzv

#### DENOISING

# Denoise the demultiplexed sequences with DADA2 
# Specify the truncation length
qiime dada2 denoise-single \
  --i-demultiplexed-seqs ./demux_seqs.qza \
  --p-trim-left 0 \
  --p-trunc-len 220 \
  --o-representative-sequences ./rep-seqs_220.qza \
  --o-table ./table_220.qza \
  --o-denoising-stats ./stats_220.qza
  
# Tabulate the denoising statistics to generate a visualization file
qiime metadata tabulate \
  --m-input-file ./stats_220.qza \
  --o-visualization ./stats_220.qzv
  
# Summarize the features table to generate a visualization file 
qiime feature-table summarize \
  --i-table ./table_220.qza \
  --o-visualization ./table_220.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt
  
# Tabulate the representative sequences to generate a visualization file
qiime feature-table tabulate-seqs \
  --i-data ./rep-seqs_220.qza \
  --o-visualization ./rep-seqs_220.qzv
  

#### GENERATION OF PHYLOGENETIC TREES

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs_220.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza


#### TAXANOMIC ASSIGNMENT

# Import a Naive Bayes Classifier pre-trained with Silva 138 99% OTUs to recognize the 515F/806R region of 16S rRNA

# Assign taxonomy to ASV’s using the pre-trained classifier
qiime feature-classifier classify-sklearn \
  --i-reads ./rep-seqs_220.qza \
  --i-classifier ./silva-138-99-515-806-nb-classifier.qza \
  --o-classification ./taxonomy.qza

# Tabulate the taxonomy results as a visualization file 
qiime metadata tabulate \
  --m-input-file ./taxonomy.qza \
  --o-visualization ./taxonomy.qzv
  

#### FILTER THE DATASET

# Filter to remove rare ASV's accounting for less than 0.005% of total reads
qiime feature-table filter-features \
  --i-table table_220.qza \
  --p-min-frequency 590 \
  --o-filtered-table frequency-filtered-table.qza

# Generate visualization file of frequency-filtered table
qiime feature-table summarize \
  --i-table frequency-filtered-table.qza \
  --o-visualization frequency-filtered-table.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt

# Filter features table by surface samples
qiime feature-table filter-samples \
  --i-table frequency-filtered-table.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --p-where "[sample_type]= 'surface'" \
  --o-filtered-table surface-freq-filtered-table.qza

# Generate visualization file of surface filtered table
qiime feature-table summarize \
  --i-table surface-freq-filtered-table.qza \
  --o-visualization surface-freq-filtered-table.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt
  
# Filter surface filtered table by sampling days before and after 3 resupply events
qiime feature-table filter-samples \
  --i-table surface-freq-filtered-table.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --p-where "[time_sampling_days] IN ('14','28','42','56','98','115')" \
  --o-filtered-table three-resupply-surface-freq-filtered-table.qza
  
# Generate visualization file of resupply event filtered table
qiime feature-table summarize \
  --i-table three-resupply-surface-freq-filtered-table.qza \
  --o-visualization three-resupply-surface-freq-filtered-table.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt
  
# Filter surface filtered table by surface material (original environmental material - plastic or wood)
qiime feature-table filter-samples \
  --i-table surface-freq-filtered-table.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --p-where "[orig_env_material] IN ('plastic','wood')" \
  --o-filtered-table plastic-wood-freq-filtered-table.qza

# Generate visualization file of surface material filtered table
qiime feature-table summarize \
  --i-table plastic-wood-freq-filtered-table.qza \
  --o-visualization plastic-wood-freq-filtered-table.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt
  
# Filter features table by skin samples
qiime feature-table filter-samples \
  --i-table frequency-filtered-table.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --p-where "[sample_type]= 'skin'" \
  --o-filtered-table skin-freq-filtered-table.qza

# Generate visualization file of skin filtered table
qiime feature-table summarize \
  --i-table skin-freq-filtered-table.qza \
  --o-visualization skin-freq-filtered-table.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt
  
# Filter skin filtered table by shower timing (before, after, exact)
qiime feature-table filter-samples \
  --i-table skin-freq-filtered-table.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --p-where "[showers_timing] IN ('before', 'after', 'exact')" \
  --o-filtered-table shower-timing-freq-filtered-table.qza

# Generate visualization file of shower timing filtered table
qiime feature-table summarize \
  --i-table shower-timing-freq-filtered-table.qza \
  --o-visualization shower-timing-freq-filtered-table.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt


#### FILTER THE DATASET PART 2

# Remove mitochondrial sequences from the resupply event filtered table
qiime taxa filter-table \
--i-table ./three-resupply-surface-freq-filtered-table.qza \
--i-taxonomy ./taxonomy.qza \
--p-exclude mitochondria \
--o-filtered-table ./resupply-table-no-mito.qza
 
# Generate visualization file of the resupply event filtered table with mitochondrial sequences removed
qiime feature-table summarize \
  --i-table resupply-table-no-mito.qza \
  --o-visualization resupply-table-no-mito.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt

# Remove mitochondrial sequences from the surface material filtered table
qiime taxa filter-table \
--i-table ./plastic-wood-freq-filtered-table.qza \
--i-taxonomy ./taxonomy.qza \
--p-exclude mitochondria \
--o-filtered-table ./plastic-wood-table-no-mito.qza

# Generate visualization file of the surface material filtered table with mitochondrial sequences removed
qiime feature-table summarize \
  --i-table plastic-wood-table-no-mito.qza \
  --o-visualization plastic-wood-table-no-mito.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt
  
# Remove mitochondrial sequences from the shower timing filtered table
qiime taxa filter-table \
--i-table ./shower-timing-freq-filtered-table.qza \
--i-taxonomy ./taxonomy.qza \
--p-exclude mitochondria \
--o-filtered-table ./shower-table-no-mito.qza

# Generate visualization file of the shower timing filtered table with mitochondrial sequences removed
qiime feature-table summarize \
  --i-table shower-table-no-mito.qza \
  --o-visualization shower-table-no-mito.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt
  
# Filter shower timing table with mitochondrial sequences removed by crew members 32 and 34
qiime feature-table filter-samples \
  --i-table shower-table-no-mito.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --p-where "[crew_id] IN ('32')" \
  --o-filtered-table shower-table-no-mito-c32.qza

qiime feature-table filter-samples \
  --i-table shower-table-no-mito.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --p-where "[crew_id] IN ('34')" \
  --o-filtered-table shower-table-no-mito-c34.qza

# Generate visualization files of the shower timing table filtered on crew members 32 and 34 
qiime feature-table summarize \
  --i-table shower-table-no-mito-c32.qza \
  --o-visualization shower-table-no-mito-c32.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt

qiime feature-table summarize \
  --i-table shower-table-no-mito-c34.qza \
  --o-visualization shower-table-no-mito-c34.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt


#### ALPHA RAREFACTION WITH FILTERED DATASETS 

# Note: for generation of alpha rarefaction curves with the resupply event-filtered features table, a new metadata file was used in 
#   which the time_sampling_days category was changed from a numerical to a categorical variable

# Alpha rarefaction of resupply event data with mitochondrial sequences removed
qiime diversity alpha-rarefaction \
  --i-table resupply-table-no-mito.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 92453 \
  --m-metadata-file hiseas_metadata2.txt \
  --o-visualization alpha-rarefaction-resupply-no-mito.qzv

# Alpha rarefaction of surface material data with mitochondrial sequences removed
qiime diversity alpha-rarefaction \
  --i-table plastic-wood-table-no-mito.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 92507 \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --o-visualization alpha-rarefaction-plastic-wood-no-mito.qzv
  
# Alpha rarefaction of shower timing data with mitochondrial sequences removed
qiime diversity alpha-rarefaction \
  --i-table shower-table-no-mito.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 120271 \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --o-visualization alpha-rarefaction-shower-no-mito.qzv
  
# Input the alpha rarefaction visualization files on q2view (https://view.qiime2.org/) to determine an appropriate rarefaction depth for each filtered dataset
  # Rarefaction depth for the resuply event filtered table with mitochondrial sequences removed - 20706 reads/sample
  # Rarefaction depth for the surface material filtered table with mitochondrial sequences removed - 30000 reads/sample
  # Rarefaction depth for the shower timing filtered table with mitochondrial sequences removed - 40000 reads/sample
  

#### COMPUTATION AND STATISTICAL ANALYSIS OF DIVERSITY METRICS

# Note: for all diversity analyses pertaining to resupply events, a new metadata file was used in 
#   which the time_sampling_days category was changed from a numerical to a categorical variable

# Calculate alpha and beta diversity metrics for the features table filtered on resupply events with mitochondrial sequences removed
# Specify a sampling depth of 20706 reads/sample 
qiime diversity core-metrics-phylogenetic \
  --i-table resupply-table-no-mito.qza \
  --i-phylogeny rooted-tree.qza \
  --m-metadata-file hiseas_metadata2.txt \
  --p-sampling-depth 20706 \
  --output-dir resupply-no-mito-core-metrics-results

# Calculate alpha group significance with each alpha diversity metric computed for resupply event data with mitochondrial sequences removed

# Faith's phylogenetic distance
qiime diversity alpha-group-significance \
  --i-alpha-diversity resupply-no-mito-core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file hiseas_metadata2.txt \
  --o-visualization resupply-no-mito-core-metrics-results/resupply-no-mito-faiths_pd_statistics.qzv

# Peilou's evenness
qiime diversity alpha-group-significance \
  --i-alpha-diversity resupply-no-mito-core-metrics-results/evenness_vector.qza \
  --m-metadata-file hiseas_metadata2.txt \
  --o-visualization resupply-no-mito-core-metrics-results/resupply-no-mito-evenness_statistics.qzv

# Observed features
qiime diversity alpha-group-significance \
  --i-alpha-diversity resupply-no-mito-core-metrics-results/observed_features_vector.qza \
  --m-metadata-file hiseas_metadata2.txt \
  --o-visualization resupply-no-mito-core-metrics-results/resupply-no-mito-observed_statistics.qzv

# Shannon diversity 
qiime diversity alpha-group-significance \
  --i-alpha-diversity resupply-no-mito-core-metrics-results/shannon_vector.qza \
  --m-metadata-file hiseas_metadata2.txt \
  --o-visualization resupply-no-mito-core-metrics-results/resupply-no-mito-shannon_statistics.qzv

# Calculate beta group significance with each beta diversity metric computed for resupply event data with mitochondrial sequences removed
# Specify the time_sampling_days metadata column to compare beta diversity of surfaces samples obtained before and after each resupply event 

# Unweighted UniFrac distance
qiime diversity beta-group-significance \
  --i-distance-matrix resupply-no-mito-core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file hiseas_metadata2.txt \
  --m-metadata-column time_sampling_days \
  --o-visualization resupply-no-mito-core-metrics-results/resupply-no-mito-unweighted-unifrac-significance.qzv \
  --p-pairwise

# Weighted UniFrac distance
qiime diversity beta-group-significance \
  --i-distance-matrix resupply-no-mito-core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file hiseas_metadata2.txt \
  --m-metadata-column time_sampling_days \
  --o-visualization resupply-no-mito-core-metrics-results/resupply-no-mito-weighted-unifrac-significance.qzv \
  --p-pairwise

# Jaccard distance
qiime diversity beta-group-significance \
  --i-distance-matrix resupply-no-mito-core-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file hiseas_metadata2.txt \
  --m-metadata-column time_sampling_days \
  --o-visualization resupply-no-mito-core-metrics-results/resupply-no-mito-jaccard-significance.qzv \
  --p-pairwise

# Bray-Curtis distance
qiime diversity beta-group-significance \
  --i-distance-matrix resupply-no-mito-core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file hiseas_metadata2.txt \
  --m-metadata-column time_sampling_days \
  --o-visualization resupply-no-mito-core-metrics-results/resupply-no-mito-bray-significance.qzv \
  --p-pairwise

# Calculate alpha and beta diversity metrics for the features table filtered on surface materials with mitochondrial sequences removed
# Specify a sampling depth of 30000 reads/sample 
qiime diversity core-metrics-phylogenetic \
  --i-table plastic-wood-table-no-mito.qza \
  --i-phylogeny rooted-tree.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --p-sampling-depth 30000\
  --output-dir material-no-mito-core-metrics-results

# Calculate alpha group significance with each alpha diversity metric computed for surface material data with mitochondrial sequences removed

# Faith's phylogenetic distance
qiime diversity alpha-group-significance \
  --i-alpha-diversity material-no-mito-core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --o-visualization material-no-mito-core-metrics-results/material-no-mito-faiths_pd_statistics.qzv

# Peilou's evenness
qiime diversity alpha-group-significance \
  --i-alpha-diversity material-no-mito-core-metrics-results/evenness_vector.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --o-visualization material-no-mito-core-metrics-results/material-no-mito-evenness_statistics.qzv

# Observed features
qiime diversity alpha-group-significance \
  --i-alpha-diversity material-no-mito-core-metrics-results/observed_features_vector.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --o-visualization material-no-mito-core-metrics-results/material-no-mito-observed_statistics.qzv

# Shannon diversity 
qiime diversity alpha-group-significance \
  --i-alpha-diversity material-no-mito-core-metrics-results/shannon_vector.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --o-visualization material-no-mito-core-metrics-results/material-no-mito-shannon_statistics.qzv

# Calculate beta group significance with each beta diversity metric computed for surface material data with mitochondrial sequences removed
# Specify the orig_env_material metadata column to compare beta diversity of surfaces samples obtained from plastic and wood surfaces

# Unweighted Unifrac
qiime diversity beta-group-significance \
  --i-distance-matrix material-no-mito-core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --m-metadata-column orig_env_material \
  --o-visualization material-no-mito-core-metrics-results/material-no-mito-unweighted-unifrac-significance.qzv \
  --p-pairwise

# Weighted Unifrac
qiime diversity beta-group-significance \
  --i-distance-matrix material-no-mito-core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --m-metadata-column orig_env_material \
  --o-visualization material-no-mito-core-metrics-results/material-no-mito-weighted-unifrac-significance.qzv \
  --p-pairwise

# Jaccard distance
qiime diversity beta-group-significance \
  --i-distance-matrix material-no-mito-core-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --m-metadata-column orig_env_material \
  --o-visualization material-no-mito-core-metrics-results/material-no-mito-jaccard-significance.qzv \
  --p-pairwise

# Bray-Curtis distance
qiime diversity beta-group-significance \
  --i-distance-matrix material-no-mito-core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --m-metadata-column orig_env_material \
  --o-visualization material-no-mito-core-metrics-results/material-no-mito-bray-significance.qzv \
  --p-pairwise

# Calculate alpha and beta diversity metrics for the shower timing features tables filtered on crew members 32 and 34
# Specify a sampling depth of 40000 reads/sample 
qiime diversity core-metrics-phylogenetic \
  --i-table shower-table-no-mito-c32.qza \
  --i-phylogeny rooted-tree.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --p-sampling-depth 40000\
  --output-dir shower-no-mito-core-metrics-results-c32

qiime diversity core-metrics-phylogenetic \
  --i-table shower-table-no-mito-c34.qza \
  --i-phylogeny rooted-tree.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --p-sampling-depth 40000\
  --output-dir shower-no-mito-core-metrics-results-c34

# Calculate alpha group significance with each alpha diversity metric computed for the shower timing data for crew members 32 and 34 with mitochondrial sequences removed

# Faith's phylogenetic distance - crew member 32
qiime diversity alpha-group-significance \
  --i-alpha-diversity shower-no-mito-core-metrics-results-c32/faith_pd_vector.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --o-visualization shower-no-mito-core-metrics-results-c32/shower-no-mito-faiths_pd_statistics-c32.qzv

# Faith's phylogenetic distance - crew member 34
qiime diversity alpha-group-significance \
  --i-alpha-diversity shower-no-mito-core-metrics-results-c34/faith_pd_vector.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --o-visualization shower-no-mito-core-metrics-results-c34/shower-no-mito-faiths_pd_statistics-c34.qzv

# Peilou's evenness - crew member 32
qiime diversity alpha-group-significance \
  --i-alpha-diversity shower-no-mito-core-metrics-results-c32/evenness_vector.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --o-visualization shower-no-mito-core-metrics-results-c32/shower-no-mito-evenness_statistics-c32.qzv

# Peilou's evenness - crew member 34
qiime diversity alpha-group-significance \
  --i-alpha-diversity shower-no-mito-core-metrics-results-c34/evenness_vector.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --o-visualization shower-no-mito-core-metrics-results-c34/shower-no-mito-evenness_statistics-c34.qzv

# Observed features - crew member 32
qiime diversity alpha-group-significance \
  --i-alpha-diversity shower-no-mito-core-metrics-results-c32/observed_features_vector.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --o-visualization shower-no-mito-core-metrics-results-c32/shower-no-mito-observed_statistics-c32.qzv

# Observed features - crew member 34
qiime diversity alpha-group-significance \
  --i-alpha-diversity shower-no-mito-core-metrics-results-c34/observed_features_vector.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --o-visualization shower-no-mito-core-metrics-results-c34/shower-no-mito-observed_statistics-c34.qzv

# Shannon diversity - crew member 32
qiime diversity alpha-group-significance \
  --i-alpha-diversity shower-no-mito-core-metrics-results-c32/shannon_vector.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --o-visualization shower-no-mito-core-metrics-results-c32/shower-no-mito-shannon_statistics-c32.qzv

# Shannon diversity - crew member 34
qiime diversity alpha-group-significance \
  --i-alpha-diversity shower-no-mito-core-metrics-results-c34/shannon_vector.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --o-visualization shower-no-mito-core-metrics-results-c34/shower-no-mito-shannon_statistics-c34.qzv

# Calculate beta group significance with each beta diversity metric computed for the shower timing data for crew members 32 and 34 with mitochondrial sequences removed
# Specify the showers_timing metadata column to compare beta diversity of skin samples obtained before and after showering

# Unweighted UniFrac - crew member 32
qiime diversity beta-group-significance \
  --i-distance-matrix shower-no-mito-core-metrics-results-c32/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --m-metadata-column showers_timing \
  --o-visualization shower-no-mito-core-metrics-results-c32/shower-no-mito-unweighted-unifrac-significance-c32.qzv \
  --p-pairwise

# Unweighted UniFrac - crew member 34
qiime diversity beta-group-significance \
  --i-distance-matrix shower-no-mito-core-metrics-results-c34/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --m-metadata-column showers_timing \
  --o-visualization shower-no-mito-core-metrics-results-c34/shower-no-mito-unweighted-unifrac-significance-c34.qzv \
  --p-pairwise

# Weighted UniFrac - crew member 32
qiime diversity beta-group-significance \
  --i-distance-matrix shower-no-mito-core-metrics-results-c32/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --m-metadata-column showers_timing \
  --o-visualization shower-no-mito-core-metrics-results-c32/shower-no-mito-weighted-unifrac-significance-c32.qzv \
  --p-pairwise

# Weighted UniFrac - crew member 34
qiime diversity beta-group-significance \
  --i-distance-matrix shower-no-mito-core-metrics-results-c34/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --m-metadata-column showers_timing \
  --o-visualization shower-no-mito-core-metrics-results-c34/shower-no-mito-weighted-unifrac-significance-c34.qzv \
  --p-pairwise

# Jaccard distance - crew member 32
qiime diversity beta-group-significance \
  --i-distance-matrix shower-no-mito-core-metrics-results-c32/jaccard_distance_matrix.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --m-metadata-column showers_timing \
  --o-visualization shower-no-mito-core-metrics-results-c32/shower-no-mito-jaccard-significance-c32.qzv \
  --p-pairwise

# Jaccard distance - crew member 34
qiime diversity beta-group-significance \
  --i-distance-matrix shower-no-mito-core-metrics-results-c34/jaccard_distance_matrix.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --m-metadata-column showers_timing \
  --o-visualization shower-no-mito-core-metrics-results-c34/shower-no-mito-jaccard-significance-c34.qzv \
  --p-pairwise

# Bray-Curtis - crew member 32
qiime diversity beta-group-significance \
  --i-distance-matrix shower-no-mito-core-metrics-results-c32/bray_curtis_distance_matrix.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --m-metadata-column showers_timing \
  --o-visualization shower-no-mito-core-metrics-results-c32/shower-no-mito-bray-significance-c32.qzv \
  --p-pairwise

# Bray-Curtis - crew member 34
qiime diversity beta-group-significance \
  --i-distance-matrix shower-no-mito-core-metrics-results-c34/bray_curtis_distance_matrix.qza \
  --m-metadata-file /mnt/datasets/project_2/hiseas/hiseas_metadata.txt \
  --m-metadata-column showers_timing \
  --o-visualization shower-no-mito-core-metrics-results-c34/shower-no-mito-bray-significance-c34.qzv \
  --p-pairwise


#### EXPORT DATA FROM QIIME2 FOR FURTHER ANALYSIS ON R

# Export the resupply event filtered features table with mitochondrial sequences removed as a biom file
qiime tools export \
--input-path resupply-table-no-mito.qza \
--output-path resupply-exported 

# Export the surface material filtered features table with mitochondrial sequences removed as a biom file
qiime tools export \
--input-path plastic-wood-table-no-mito.qza \
--output-path plastic-wood-exported

# Export the shower timing filtered features table with mitochondrial sequences removed as a biom file
qiime tools export \
--input-path shower-table-no-mito.qza \
--output-path shower-exported 

# Export the taxonomic classification for each ASV as a tsv file
qiime tools export \
--input-path taxonomy.qza \
--output-path shower-exported 

# Export the rooted phylogenetic tree as a .nwk file 
qiime tools export \
--input-path rooted-tree.qza \
--output-path shower-exported

# Manually edit column names of exported taxonomy file using nano
nano shower-exported/taxonomy.tsv
  # 'Feature ID' to '#OTUID'
  # 'Taxon' to 'taxonomy'
  # 'Confidence' to 'confidence'

# Combine taxonomy with BIOM data

biom add-metadata \
-i resupply-exported/feature-table.biom \
-o resupply-exported/resupply-table-with-taxonomy.biom \
--observation-metadata-fp shower-exported/taxonomy.tsv \
--sc-separated taxonomy

biom add-metadata \
-i plastic-wood-exported/feature-table.biom \
-o plastic-wood-exported/plastic-wood-table-with-taxonomy.biom \
--observation-metadata-fp shower-exported/taxonomy.tsv \
--sc-separated taxonomy

biom add-metadata \
-i shower-exported/feature-table.biom \
-o shower-exported/shower-table-with-taxonomy.biom \
--observation-metadata-fp shower-exported/taxonomy.tsv \
--sc-separated taxonomy




