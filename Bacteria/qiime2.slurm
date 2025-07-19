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
  --p-trim-left-f 21 \
  --p-trim-left-r 23 \
  --p-trunc-len-f 248 \
  --p-trunc-len-r 233 \
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

# Downloading the SILVA database files
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/taxonomy/tax_slv_ssu_138.2.txt.gz
gunzip tax_slv_ssu_138.2.txt.gz

wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.2.txt.gz
gunzip taxmap_slv_ssu_ref_nr_138.2.txt.gz

wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/taxonomy/tax_slv_ssu_138.2.tre.gz
gunzip tax_slv_ssu_138.2.tre.gz

wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz
gunzip SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz

# Converting SILVA files to QIIME2 artifacts 
qiime tools import \
    --type 'FeatureData[SILVATaxonomy]' \
    --input-path tax_slv_ssu_138.1.txt \
    --output-path taxranks-silva-138.1-ssu-nr99.qza
    
qiime tools import \
    --type 'FeatureData[SILVATaxidMap]' \
    --input-path taxmap_slv_ssu_ref_nr_138.1.txt \
    --output-path taxmap-silva-138.1-ssu-nr99.qza
    
qiime tools import \
    --type 'Phylogeny[Rooted]' \
    --input-path tax_slv_ssu_138.1.tre \
    --output-path taxtree-silva-138.1-nr99.qza

qiime tools import \
    --type 'FeatureData[RNASequence]' \
    --input-path SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta \
    --output-path silva-138.1-ssu-nr99-rna-seqs.qza

# Convert SILVA sequences to DNA sequences
qiime rescript reverse-transcribe \
    --i-rna-sequences silva-138.1-ssu-nr99-rna-seqs.qza \
    --o-dna-sequences silva-138.1-ssu-nr99-seqs.qza

# Creating SILVA reference database 
qiime rescript parse-silva-taxonomy \
    --i-taxonomy-tree taxtree-silva-138.1-nr99.qza \
    --i-taxonomy-map taxmap-silva-138.1-ssu-nr99.qza \
    --i-taxonomy-ranks taxranks-silva-138.1-ssu-nr99.qza \
    --o-taxonomy silva-138.1-ssu-nr99-tax.qza

# “Culling” low-quality sequences with cull-seqs
qiime rescript cull-seqs \
    --i-sequences silva-138.1-ssu-nr99-seqs.qza \
    --o-clean-sequences silva-138.1-ssu-nr99-seqs-cleaned.qza

# Filtering sequences by length and taxonomy
qiime rescript filter-seqs-length-by-taxon \
    --i-sequences silva-138.1-ssu-nr99-seqs-cleaned.qza \
    --i-taxonomy silva-138.1-ssu-nr99-tax.qza \
    --p-labels Archaea Bacteria Eukaryota \
    --p-min-lens 900 1200 1400 \
    --o-filtered-seqs silva-138.1-ssu-nr99-seqs-filt.qza \
    --o-discarded-seqs silva-138.1-ssu-nr99-seqs-discard.qza

# Dereplication of sequences and taxonomy
qiime rescript dereplicate \
    --i-sequences silva-138.1-ssu-nr99-seqs-filt.qza  \
    --i-taxa silva-138.1-ssu-nr99-tax.qza \
    --p-mode 'uniq' \
    --o-dereplicated-sequences silva-138.1-ssu-nr99-seqs-derep-uniq.qza \
    --o-dereplicated-taxa silva-138.1-ssu-nr99-tax-derep-uniq.qza

# Creating the classifier 
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads  silva-138.1-ssu-nr99-seqs-derep-uniq.qza \
  --i-reference-taxonomy silva-138.1-ssu-nr99-tax-derep-uniq.qza \
  --o-classifier silva-138.1-ssu-nr99-classifier.qza

# Creating an amplicon-region-specific classifier
qiime feature-classifier extract-reads \
    --i-sequences silva-138.1-ssu-nr99-seqs-derep-uniq.qza \
    --p-f-primer GACTCCTACGGGAGGCWGCAG \
    --p-r-primer GACTACHVGGGTATCTAATCC \
    --p-n-jobs 2 \
    --p-read-orientation 'forward' \
    --o-reads silva-138.1-ssu-nr99-seqs-337f-805r.qza

# Dereplicate the extracted region
qiime rescript dereplicate \
    --i-sequences silva-138.1-ssu-nr99-seqs-337f-805r.qza \
    --i-taxa silva-138.1-ssu-nr99-tax-derep-uniq.qza \
    --p-mode 'uniq' \
    --o-dereplicated-sequences silva-138.1-ssu-nr99-seqs-337f-805r-uniq.qza \
    --o-dereplicated-taxa  silva-138.1-ssu-nr99-tax-337f-805r-derep-uniq.qza

# Creating the amplicon-region specific classifier
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads silva-138.1-ssu-nr99-seqs-337f-805r-uniq.qza \
    --i-reference-taxonomy silva-138.1-ssu-nr99-tax-337f-805r-derep-uniq.qza \
    --o-classifier silva-138.1-ssu-nr99-337f-805r-classifier.qza

# Taxonomic classification of ASV
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138.1-ssu-nr99-337f-805r-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

# Convert taxonomy artifact to a visualization file.
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

#Remove mitochondria and chloroplast
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
  --p-max-depth 90000 \
  --m-metadata-file metadata.tsv \
  --o-visualization alpha-rarefaction.qzv

# Calculate diversity metrics. 
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 13700 \
  --m-metadata-file metadata-5Resistant.tsv \
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
# Filtering ASV with a relative abundance of 0.01% from the total reads
qiime feature-table filter-features \
  --i-table table.qza \
  --p-min-frequency 158 \
  --o-filtered-table filtered-table.qza

# Collapsing table to species level
echo "get taxa for filtered table" 
qiime taxa collapse \
  --i-table filtered-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table table-l7-filtered.qza

# Running the ANCOM-BC analysis
qiime composition ancombc \
  --i-table table-l7-filtered.qza \
  --p-lib-cut 13700 \
  --m-metadata-file metadata2.tsv \
  --p-formula 'Resistance' \
  --o-differentials l7-ancombc-resistance-filtered-5R.qza

# Creating the ANCOM-BC plot
qiime composition da-barplot \
  --i-data l7-ancombc-resistance-filtered-5R.qza \
  --p-significance-threshold 0.001 \
  --p-level-delimiter ';' \
  --o-visualization l7-da-barplot-resistance-filtered-5R.qzv  


echo "End time"
date +"%d %B %Y %H:%M:%S"

