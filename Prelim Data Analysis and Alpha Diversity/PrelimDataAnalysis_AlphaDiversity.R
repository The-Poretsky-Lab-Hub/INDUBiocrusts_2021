### Preliminary data analysis and alpha diversity


########################## Import Packages #####################################

library(dada2)
library(ggplot2)
library(Biostrings)
library(phyloseq)
library(plyr)
library("phylosmith")
library("DESeq2")
library(ggpubr)
library(ggeasy)
library(vegan)
library("microbiome")
library(microbiomeMarker)
library("tidyverse")
library("ape")
library("dendextend")
library(pairwiseAdonis)
library('microeco')
library("ALDEx2")
library("microViz")
library(randomForest)
library(knitr)
library(doBy)
library(plotrix)
library(reshape2)
library(doParallel)
library(ggh4x)

######################## Load Previous DADA2 results ########################

#read ASV table and taxa table
st <- readRDS("~/OneDrive/Documents_OneDrive/PHD/RESEARCH/BIOCRUSTS/biocrust_st_10May2023.rds")

taxtab <- readRDS("~/OneDrive/Documents_OneDrive//PHD/RESEARCH/BIOCRUSTS/biocrust_taxtab_10May2023.rds")

### Handoff to Phyloseq [.]
#import dataframe
#can add metadata here
samdf <- read.csv("~/OneDrive/Documents_OneDrive/PHD/RESEARCH/BIOCRUSTS/BIOCRUST SAMPLES/ALL SAMPLES/BIOCRUST SAMPLES/Biocrust_Metadata_Nutrients.csv",
                  row.names=1)

#construct phyloseq object
ps <- phyloseq(otu_table(st, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxtab))

ps
ps <- subset_taxa(ps, Kingdom =="Bacteria") 
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
ps <- subset_taxa(ps, (Class!="Chloroplast") | is.na(Class))
ps <- subset_taxa(ps, (Family!="Mitochondria"))


sample_data(ps)$Distance.from.Shore.m. <- as.factor(sample_data(ps)$Distance.from.Shore.m.)

sample_data(ps)$Read_Count <- sample_sums(ps)

sample_data(ps)



head(taxa_names(ps))
head(taxa_sums(ps))
mean(sample_sums(ps))

# calculate mean species richness per sample
ps.alpha.df <- as.data.frame(richness(ps))
mean(ps.alpha.df$observed)

########################## Reorder categories ##################################
sample_names(ps) <- sample_data(ps)$ID

sample_data(ps)$Dune.Morphology <- factor(sample_data(ps)$Dune.Morphology,     # Reorder factor levels
                                          c("Linear", "Parabolic", "Oak Savannah"))
sample_data(ps)$Sample.Type <- factor(sample_data(ps)$Sample.Type,     # Reorder factor levels
                                      c("Bare Soil", "Sub-Crust", "Biocrust"))
ps <- ps_arrange(ps, target = sample_data(ps)$Sample.Type)
sample_data(ps)

order <- rownames(sample_data(ps))
ps <- ps_reorder(ps, sample_order = order)

taxa_names(ps) <- paste0("Seq", seq(ntaxa(ps)))


################################ Remove NA Samples #############################
ps <- subset_samples(ps, sample_data(ps)$Dune.Morphology != 'Oak Savannah')
ps <- subset_samples(ps, sample_data(ps)[,6:12] != 'NA')
ps <- subset_samples(ps, sample_data(ps)$ID != '7CT')
ps <- subset_samples(ps, sample_data(ps)$ID != '10C-T1')
ps <- subset_samples(ps, sample_data(ps)$ID != '11A1')
ps <- subset_samples(ps, sample_data(ps)$ID != '11A2')
ps <- subset_samples(ps, sample_data(ps)$ID != '12B1')
ps
sample_data(ps)

######################### Data transformations #################################
#Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))

# CLR Transformation 
ps.clr <- 
  transform(
    ps,
    transform = "clr",
    target = "OTU",
    shift = 0,
    scale = 1,
    log10 = FALSE,
    reference = 1)

######################## Preliminary Data Analysis with Phyloseq ###############

###                       Define Comparisons & Colors
my_comparisons <- list( c("Bare Soil", "Biocrust"), c("Bare Soil", "Sub-Crust"), c("Sub-Crust", "Biocrust"))

colors <- c('darkgoldenrod1', 'cadetblue3', 'lightgreen', 'darkgoldenrod1', 'cadetblue3', 'lightgreen')

###Richness plots
###                           Sample Type
p1 <-
  plot_richness(ps, x="Sample.Type", measures=c("Observed", "Shannon"))+
  geom_boxplot(alpha = 0.5, width = 0.5, 
               color = "black", fill = colors)+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))+
  theme_bw()+
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(method = "kruskal.test", label.x = 1.5, size = 4, label.y = 0)+
  rotate_x_text(angle = 360, hjust = 0.5, size = 12)+
  ylab("Alpha Diversity")+
  xlab("Sample Type")+
  theme(text=element_text(size=13))+
  ylim(limits = c(0,NA))


###                     Average Read Count
###                       Define Comparisons & Colors
my_comparisons <- list( c("Bare Soil", "Biocrust"), c("Bare Soil", "Sub-Crust"), c("Sub-Crust", "Biocrust"))
colors <- c('darkgoldenrod1', 'cadetblue3', 'lightgreen')

p2 <- 
  ggplot(data.frame(sample_data(ps)), aes(x = Sample.Type, y = Read_Count))+
  geom_boxplot(col = "black", fill = colors, alpha = 0.5)+
  geom_point()+
  theme_pubr()+
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(method = "kruskal.test", label.x = 1.7, label.y = 0, size = 4.5)+
  ylab("Number of Merged Reads")+
  xlab("Sample Type")+
  theme(text=element_text(size=13))+
  ylim(limits = c(0,NA))


prelim <- ggarrange(p1, p2, labels = "AUTO")

#Save Plots 
ggsave("~/OneDrive/Documents_OneDrive/PHD/RESEARCH/DADA2/Biocrust Figures/For_Publication/ prelim_reduced_updated_20Aug.png",
       prelim, width = 12, height = 6)

