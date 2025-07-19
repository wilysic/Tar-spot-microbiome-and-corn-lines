#!/bin/bash
#SBATCH -A highmem
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --job-name=dada-taxonomy
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --error=%x-%J-%u.err
#SBATCH --output=%x-%J-%u.out

echo "Start time"
date +"%d %B %Y %H:%M:%S"

# Loading QIIME2
echo "loading anaconda"
module load anaconda/2022.10-py39

echo "loading qiime2-2023.9"
conda activate qiime2-amplicon-2024.2

# Import sequencing data
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path pe-33-manifest.tsv \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

# Visualizing imported data
qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization paired-end-demux.qzv

# Quality control and clustering sequences in ASV
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux.qza \
  --p-trim-left-f 22 \
  --p-trim-left-r 21 \
  --p-trunc-len-f 190 \
  --p-trunc-len-r 190 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza

#Convert .qza to .qzv format - visualization files
qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

# Creating the Fungi database
# Location of UNITED developer files: https://dx.doi.org/10.15156/BIO/2938080
qiime tools import \
--type FeatureData[Sequence] \
--input-path sh_refs_qiime_ver9_99_s_25.07.2023_dev.fasta \
--output-path united_refs_ver9_99_s_25.07.2023_dev.qza

qiime tools import \
--type FeatureData[Taxonomy] \
--input-path sh_taxonomy_qiime_ver9_99_s_25.07.2023_dev.txt \
--output-path united_taxonomy_ver9_99_s_25.07.2023_dev.qza \
--input-format HeaderlessTSVTaxonomyFormat

qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads united_refs_ver9_99_s_25.07.2023_dev.qza \
--i-reference-taxonomy united_taxonomy_ver9_99_s_25.07.2023_dev.qza \
--o-classifier unite-ver9-wsf-99-classifier-25.07.2023.qza

# Taxonomic classification of ASV
qiime feature-classifier classify-sklearn \
  --i-classifier unite-ver9-wsf-99-classifier-25.07.2023.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

# Convert taxonomy artifact to a visualization file.
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

# Remove mitochondria and chloroplast 
qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table table-no-mitochondria-no-chloroplast.qza

#Rename files so that the filtered one is table.qza
mv table.qza table-with-mitochondria-and-chloroplast.qza
mv table-no-mitochondria-no-chloroplast.qza table.qza

# Creating taxonomic bar plot
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization taxa-bar-plots.qzv


# Create a phylogenetic tree.
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

# Creating a figure with rarefaction curves 
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 95000 \
  --m-metadata-file metadata.tsv \
  --o-visualization alpha-rarefaction.qzv

#Calculate diversity metrics. 
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 67200 \
  --m-metadata-file metadata.tsv \
  --output-dir core-metrics-results

# Statistical testing for alpha metrics
# Faith
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

# Evenness
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

# Shannon
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization core-metrics-results/shannon-group-significance.qzv


# Statistical testing for beta metrics
#Unweighted unifract for resistance
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column Resistance \
  --o-visualization core-metrics-results/unweighted-unifrac-resistance-significance.qzv \
  --p-pairwise

#  Bray-Curtis for resistance  
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column Resistance \
  --o-visualization core-metrics-results/bray-curtis-resistance-significance.qzv \
  --p-pairwise

# ANCOM-BC between resistant and susceptible corn lines
# Collapsing table to species level
qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table table-l7.qza

# Running the ANCOM-BC analysis
qiime composition ancombc \
  --i-table table-l7.qza \
  --p-lib-cut 60000 \
  --m-metadata-file metadata.tsv \
  --p-formula 'Resistance' \
  --o-differentials l7-ancombc-resistance-4R.qza

# Creating the ANCOM-BC plot
qiime composition da-barplot \
  --i-data l7-ancombc-resistance-4R.qza \
  --p-significance-threshold 0.001 \
  --p-level-delimiter ';' \
  --o-visualization l7-da-barplot-resistance-4R.qzv

echo "End time"
date +"%d %B %Y %H:%M:%S"

