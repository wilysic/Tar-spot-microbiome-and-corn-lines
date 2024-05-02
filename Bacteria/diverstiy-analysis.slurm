#!/bin/bash
#SBATCH -A highmem
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --job-name=diversity-analysis
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --error=%x-%J-%u.err
#SBATCH --output=%x-%J-%u.out

echo "Start time"
date +"%d %B %Y %H:%M:%S"

#step 1, cd into the proper directory. This directory must already exist

cd $RCAC_SCRATCH
cd Bacteria_RS

pwd

#step 2, Load Modules

echo "loading anaconda"
module load anaconda

echo "loading qiime2-2023.9"
conda activate qiime2-amplicon-2024.2

echo "initiating qiime2 pipeline"

#step 3, Run the code you want
echo "Core metrics"

#Calculate diversity metrics. You will need to designate a --p-sampling-depth. You will want to consult the rarifaction output to help you make a decision.
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 11600 \
  --m-metadata-file metadata.tsv \
  --output-dir core-metrics-results

qiime diversity alpha \
  --i-table table.qza \
  --p-metric chao1 \
  --o-alpha-diversity core-metrics-results/chao1_vector.qza

echo "qiime diversity alpha-group-significance"

#Statistical testing
#Faith
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

#Evenness
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

#Shannon
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization core-metrics-results/shannon-group-significance.qzv

#Observed features
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/observed_features_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization core-metrics-results/observed_features-group-significance.qzv

#Chao1
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/chao1_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization core-metrics-results/chao1-group-significance.qzv
  
echo "qiime diversity alpha-correlation"
#########
#Faith
qiime diversity alpha-correlation \
    --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
    --m-metadata-file metadata.tsv \
    --o-visualization core-metrics-results/faith-alpha-correlation.qzv

#Eveness
qiime diversity alpha-correlation \
    --i-alpha-diversity core-metrics-results/evenness_vector.qza \
    --m-metadata-file metadata.tsv \
    --o-visualization core-metrics-results/eveness-alpha-correlation.qzv

#Shannon
qiime diversity alpha-correlation \
    --i-alpha-diversity core-metrics-results/shannon_vector.qza \
    --m-metadata-file metadata.tsv \
    --o-visualization core-metrics-results/shannon-alpha-correlation.qzv

#Observed features
qiime diversity alpha-correlation \
    --i-alpha-diversity core-metrics-results/observed_features_vector.qza \
    --m-metadata-file metadata.tsv \
    --o-visualization core-metrics-results/observed_features-alpha-correlation.qzv
    
echo "qiime diversity beta group significance"
#For resistance

#Unweighted unifract for resistance
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column Resistance \
  --o-visualization core-metrics-results/unweighted-unifrac-resistance-significance.qzv \
  --p-pairwise

#weighted unifract for resistence
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column Resistance \
  --o-visualization core-metrics-results/weighted-unifrac-resistance-significance.qzv \
  --p-pairwise
  
#jaccard for resistence
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column Resistance \
  --o-visualization core-metrics-results/jaccard-resistance-significance.qzv \
  --p-pairwise

#bray-curtis for resistence  
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column Resistance \
  --o-visualization core-metrics-results/bray-curtis-resistance-significance.qzv \
  --p-pairwise


echo "End time"
date +"%d %B %Y %H:%M:%S"

