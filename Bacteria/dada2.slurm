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
echo "dada2"

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux.qza \
  --p-trim-left-f 21 \
  --p-trim-left-r 23 \
  --p-trunc-len-f 248 \
  --p-trunc-len-r 233 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza

#Generate Taxonomy.qza

echo "Classify and get .qzv"
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138.1-ssu-nr99-337f-805r-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

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

echo "Make a stacked barplot"
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization taxa-bar-plots.qzv

#Convert .qza to .qzv format

qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv


#Create a .qzv from the output of denoise
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv


echo "Create a phylogenetic tree. (Phylogenetic method)"
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

echo "Rarefaction"
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 90000 \
  --m-metadata-file metadata.tsv \
  --o-visualization alpha-rarefaction.qzv


echo "End time"
date +"%d %B %Y %H:%M:%S"

