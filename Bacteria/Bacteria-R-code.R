## This script goes through making a taxa bar plot, alpha and beta diversity plots
## Code adapted from https://github.com/john2929/ANSC516/tree/master/Tutorials

# Install phyloseq:
 if (!requireNamespace("BiocManager", quietly = TRUE))
 install.packages("BiocManager")
BiocManager::install("phyloseq")

install.packages("ggpubr")
install.packages("writexl")
install.packages("FSA")
install.packages("multcompView")
install.packages("tidyverse")
install.packages("vegan")
install.packages("devtools")
devtools::install_github("jbisanz/qiime2R")

library(qiime2R)
library(phyloseq)
library(tidyverse)
library(multcompView)
library(FSA)
library(dplyr)
library(ggpubr)
library(writexl)
library(ggplot2)
library(dunn.test)
library(devtools)
library(vegan)


getwd()
setwd("/Users/Fungi_Results")
list.files()

# Creating a folder for the outputs
if(!dir.exists("output/taxa"))
  dir.create("output/taxa")

# Creating a phyloseq object
physeq <- qza_to_phyloseq(
  features="core-metrics-results/rarefied_table.qza",
  tree="rooted-tree.qza",
  taxonomy = "taxonomy.qza",
  metadata = "metadata.tsv"
)

asv_table <- data.frame(otu_table(physeq), check.names = F)
metadata <- data.frame(sample_data(physeq), check.names = F)
taxonomy <- data.frame(tax_table(physeq), check.names = F)

levels(metadata$Corn.line)

# Clean up taxonomy
head(taxonomy)
tax.clean <- taxonomy

# ASVs unclassified at any level are classified as the lowest taxonomic available
tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("uncl_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("uncl_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("uncl_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("uncl_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("uncl_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("uncl_",tax.clean$Genus[i], sep = "_")
  }
}


# Assign tables as variables into phyloseq
OTU.physeq = otu_table(as.matrix(asv_table), taxa_are_rows=TRUE)
tax.physeq = tax_table(as.matrix(tax.clean))    
meta.physeq = sample_data(metadata)

# Merge variagles into an object of class phyloseq.
physeq_bar_plot = phyloseq(OTU.physeq, tax.physeq, meta.physeq)


# Set colors for plotting
my_colors <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#ffff99','#ff7f00','#cab2d6','#6a3d9a','#fdbf6f','#b15928', 
  "#CBD588","#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
)

# Defining taxonomic levl
my_level <- c("Phylum", "Family", "Genus")

# Plot taxa based on corn lines
my_column <- "Corn.line"  
my_column_ordered <- c("CML103", "CML69",  "685806", "685831", "685915", "685790", 
                       "685919", "685950", "685836", "685918","B97", "4401350", 
                       "685788", "685920")

# Set abundance of taxa in at least one sample to plot only 30 most abundant
abund_filter <- 0.0060  

for(ml in my_level){
  print(ml)
  taxa.summary <- physeq_bar_plot %>%
    tax_glom(taxrank = ml, NArm = FALSE) %>%  # agglomerate at `ml` level
    transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
    psmelt()  %>%                               # Melt to long format
    group_by(get(my_column), get(ml)) %>%
    summarise(Abundance.average=mean(Abundance)) 
  taxa.summary <- as.data.frame(taxa.summary)
  colnames(taxa.summary)[1] <- my_column
  colnames(taxa.summary)[2] <- ml
  physeq.taxa.max <- taxa.summary %>% 
    group_by(get(ml)) %>%
    summarise(overall.max=max(Abundance.average))
  physeq.taxa.max <- as.data.frame(physeq.taxa.max)
  colnames(physeq.taxa.max)[1] <- ml
  physeq.taxa.mean <- taxa.summary %>% 
    group_by(get(ml)) %>%
    summarise(overall.mean=mean(Abundance.average))
  physeq.taxa.mean <- as.data.frame(physeq.taxa.mean)
  colnames(physeq.taxa.mean)[1] <- ml
  # merging the phyla means with the metadata #
  physeq_meta <- merge(taxa.summary, physeq.taxa.max)
  physeq_meta <- merge (physeq_meta, physeq.taxa.mean)
  physeq_meta_filtered <- filter(physeq_meta, overall.max>abund_filter)
  #str(physeq_meta_filtered)
  physeq_meta_filtered$my_column_ordered = factor(physeq_meta_filtered[[my_column]], my_column_ordered)
  physeq_meta_filtered[[ml]] <- factor(physeq_meta_filtered[[ml]])
  y = tapply(physeq_meta_filtered$overall.mean, physeq_meta_filtered[[ml]], function(y) max(y))
  y = sort(y, TRUE)
  physeq_meta_filtered[[ml]] = factor(as.character(physeq_meta_filtered[[ml]]), levels=names(y))
  levels(physeq_meta_filtered[[ml]])
  # Creating the Taxonomic bar plot 
  ggplot(physeq_meta_filtered, aes(x = my_column_ordered, y = Abundance.average, fill = get(ml))) + 
    #facet_grid(.~.) +
    geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = my_colors) +
    # Remove x axis title
    #theme(axis.title.x = element_blank()) + 
    ylim(c(0,1)) +
    guides(fill = guide_legend(reverse = TRUE, keywidth = .5, keyheight = .5, ncol = 1)) +
    theme(legend.text=element_text(size=8)) +
    #theme(legend.position="bottom") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(legend.title = element_blank()) +
    ylab("Relative Abundance") +
    xlab(my_column) +
    ggtitle(paste0(ml, " (>", abund_filter * 100,"%) in at least 1 sample")) 
  ggsave(paste0("output/taxa/", ml, "BarPlot_", my_column, ".png"), height = 5, width = 8)
}



## ALPHA figures

meta<-read_q2metadata("metadata.tsv")
colnames(meta)[4] <- "Corn.line"

# Completing sample name 
meta$Corn.line.ord = factor(meta$Corn.line, 
                            levels = c("CML103", "CML69", "685806","685831","CML52", "685915", "685790", "685919", "685950", "685836", "685918","B97", "4401350", "685788", "685920","TX303" ),
                            labels =  c("CML103", "CML69", "PP685806", "PP685831", "CML52",  "PP685915", "PP685790", "PP685919", "PP685950", "PP685836", "PP685918","B97", "4401350", "PP685788", "PP685920", 'TXT303'))

levels(meta$Corn.line.ord)

evenness = read_qza("core-metrics-results/evenness_vector.qza")
evenness<-evenness$data %>% rownames_to_column("Corn.line.ord") 

shannon = read_qza("core-metrics-results/shannon_vector.qza")
shannon<-shannon$data %>% rownames_to_column("Corn.line.ord") 

faith_pd = read_qza("core-metrics-results/faith_pd_vector.qza")
faith_pd<-faith_pd$data %>% rownames_to_column("Corn.line.ord") 
faith_pd<-subset(faith_pd, select = -Corn.line.ord)
colnames(faith_pd)[1] <- "Corn.line.ord"
colnames(faith_pd)[2] <- "faith_pd"


# Merge alpha tables to metadata
alpha_diversity = merge(x=faith_pd, y=evenness, by.x = "Corn.line.ord", by.y = "Corn.line.ord")
alpha_diversity = merge(alpha_diversity, observed_features, by.x = "Corn.line.ord", by.y = "Corn.line.ord")
alpha_diversity = merge(alpha_diversity, shannon, by.x = "Corn.line.ord", by.y = "Corn.line.ord")
meta = merge(meta, alpha_diversity, by.x = "SampleID", by.y = "Corn.line.ord")
row.names(meta) <- meta$SampleID
str(meta)

my_colors <- c(
  '#a6cee3', '#1f78b4', '#b2df8a',  '#33a02c', '#e31a1c', '#fdbf6f', 
  '#ff7f00', '#cab2d6',  '#6a3d9a', '#ffff99', '#b15928', "#CBD588",
  "#5F7FC7", "#D14285",  "#C84248", "#8569D5", "#673770", "#652926",
  "#508578", "orange",   "#DA5724", '#fb9a99', "#CD9BCD")
  

##########
#This following part is to add the significant letter above alpha bar plots
#The code create Shannon and Faith's plot replace whit the corresponding name 

# Run Dunn test
pw <- dunnTest(meta$pielou_evenness ~ meta$Corn.line.ord, method = "bh")
print(pw)

# Extract adjusted p-values and comparison names
adjusted_pvalues <- pw$res$P.adj
comparisons <- pw$res$Comparison

# Remove spaces from comparison names
comparisons_clean <- gsub(" ", "", comparisons)

# Create a logical vector indicating significant differences
Diff <- adjusted_pvalues < 0.05
names(Diff) <- comparisons_clean

# Apply multcompLetters to obtain letter groupings
letters <- multcompLetters(Diff)$Letters

# Convert to a data frame
letter_df <- data.frame(Corn.line.ord = names(letters), LetterGroup = letters)

# Compute median values for each group (to position letters)
meta_summary <- meta %>%
  group_by(Corn.line.ord) %>%
  summarize(median_evenness = median(pielou_evenness, na.rm = TRUE)) %>%
  left_join(letter_df, by = "Corn.line.ord")


# Create the boxplot with letter annotations
evenness_boxplot <- ggplot(meta, aes(x = Corn.line.ord, y = pielou_evenness, fill = factor(Corn.line.ord))) + 
  geom_boxplot() + 
  scale_fill_manual(values = my_colors) + 
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y = "Evenness Fungi", x = "") +
  # Add letters only once per group, positioned at the median of each group
  geom_text(data = meta_summary, aes(x = Corn.line.ord, y = median_evenness, label = LetterGroup), 
            position = position_dodge(width = 0.75), vjust = -3, size = 3)
ggsave("output/Evenness_fungi.png", evenness_boxplot, height = 5, width = 5)



## BETA plots
metadata<-read_q2metadata("metadata.tsv")
colnames(metadata)[3] <- "Resistance"

# Importing data
bc_PCoA<-read_qza("core-metrics-results/bray_curtis_pcoa_results.qza")

# Assigning colors
resistance_colors <- c("Orange", "Green")

# Defining variable
my_column <- "Resistance"

# Bray Curtis
# Extracting values
bc_meta <- bc_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

# Defining centroids
centroids <- aggregate(cbind(PC1,PC2)~get(my_column),bc_meta,mean)
colnames(centroids)[1] <- "Resistance"

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=resistance_colors, name = my_column)
ggsave(paste0("output/BC-ellipse_", my_column,".png"), height=4.5, width=4.5, device="png") 


# Unweighted UniFrac
Uwuni_PCoA<-read_qza("core-metrics-results/unweighted_unifrac_pcoa_results.qza")

# Extracting values
Uwuni_meta <- Uwuni_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),Uwuni_meta,mean)
colnames(centroids)[1] <- "Resistance"

ggplot(Uwuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + 
  geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Uwuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Uwuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=resistance_colors, name = "Resistance")
ggsave(paste0("output/Uwuni-ellipse_", my_column,".png"), height=4.5, width=4.5, device="png") 

