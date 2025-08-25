######## DADA2 BIOCRUST SAMPLES ##############

# Set Up Data

setwd("~/OneDrive/Documents (OneDrive)/PHD/RESEARCH/BIOCRUSTS/BIOCRUST SAMPLES/ALL SAMPLES/BIOCRUST SAMPLES/")

.libPaths("C:/R Libraries")
.libPaths()

library(devtools)
library(dada2)
library(ggplot2)
library(Biostrings)
library(phyloseq)
library("dplyr")
library(ellipsis)

###SET PATH AND CHECK FILES###

path<- "~/OneDrive/Documents_OneDrive/PHD/RESEARCH/BIOCRUSTS/BIOCRUST SAMPLES/ALL SAMPLES/BIOCRUST SAMPLES/Manuscript Sequences/"

list.files(path)

data.frame(list.files(path))

write.csv(data.frame(list.files(path)), 
          "~/OneDrive/Documents_OneDrive/UIC Work/Manuscripts/sample_files.csv")

# 56 raw read files (1F and 1R for each sample) for 28 samples. 
# It looks like 2 SAMPLES ARE MISSING. Will have to find these another time
#Added one additional oak savannah sample and removed 5AT



###SORT INTO FORWARD AND REVERSE READS###

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
head (fnFs)
head(fnRs)


###EXTRACT SAMPLES NAMES###


sample.names <- sapply(strsplit(basename(fnFs), "_S"), `[`, 1)
head (sample.names)
sample.names
# Successully extrcted 28 sample names 

###CHECK QUALITY OF READS###

plotQualityProfile(fnFs[28])
plotQualityProfile(fnRs[28])

# Quality plots for all samples look surpisingly good. There is no need to trim these samples. 
# However, primers still need to be removed. 


# Primer: 515f (19) 806r (20) Amplicon Length: (806-515) = 291 (-(21+20) = 250)
# Updated sequences: 515F (Parada)–806R (Apprill), forward-barcoded:
# FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT
## For some reason, an additional 2 bp are occurring in front of my forward primer. Instead
# of removing 19 bp, we will remove the first 21. 

###Filter and place in sub-directory 

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

###Trim and filter reads 

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=FALSE, multithread=TRUE, trimLeft = c(21,20))

#Since both our F and R reads had score above 30 for all cycles, we will not trim our read and 
#leave them at 160 nts in length. Just to note: read length should be 20+boilogical.length.variation
# in order to be able to properly merge later on. 

##Run did not work when leaving reads untrimmed. Instead, we trimmed reads to 135 nts each. 

##Most of the reads remained in our samples, which is excellent!


###LEARN ERROR RATES 

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE) 

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

###SAMPLE INFERENCE 

dadaFs <-dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[3]]
dadaRs[[3]]


###MERGE PAIRED READS 

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE,
                      maxMismatch = 0,
                      minOverlap = 12)
                      #justConcatenate = TRUE)
head(mergers[[1]])

#For some reason, very few of my reads are merging, regardless of how I trim my reads.
#Currently, I have concatenated my reads but this may impact downstream analysis. I should be 
#critical of my interpreation of the data after this point. 

#UPDATE: No concatenation needed. Reads are able to merge when reads are no truncated. Since 
#read quality is already good, truncation really isn't necessary. 

###CONSTRUCT SEQUENCE TABLE

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# [1]     28 37918


###INSPECT DISTRIBUTION OF SEQUENCE LENGTHS

table(nchar(getSequences(seqtab)))

#REMOVE LOWER LENGTH MERGED READS 

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 248:256]

table(nchar(getSequences(seqtab2)))

###REMOVE CHIMERAS

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

#Filtered merged reads
seqtab2.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab2.nochim)

sum(seqtab2.nochim)/sum(seqtab2)
# [1] 0.9581179
# Chimeras only make up 5% of my merged reads 


###TRACK READS THROUGH THE PIPELINE 

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab2.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

track.df <- as.data.frame(track)

sum(track.df$input)
mean(track.df$nonchim)

###ASSIGN TAXONOMY 

taxa <- assignTaxonomy(seqtab2.nochim, "~/OneDrive/Documents (OneDrive)/PHD/RESEARCH/DADA2/silva_nr99_v138.1_train_set.fa.gz", 
                       multithread=TRUE,tryRC=TRUE)

#ADD SPECIES 

taxa.species <- addSpecies(taxa, "~/OneDrive/Documents (OneDrive)/PHD/RESEARCH/DADA2/silva_species_assignment_v138.1.fa.gz")

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
taxa.print


###SAVE ASV AND TAXA TABLES###

saveRDS(seqtab2.nochim, "~/OneDrive/Documents (OneDrive)/PHD/RESEARCH/BIOCRUSTS/biocrust_st_10May2023.rds")
saveRDS(taxa, "~/OneDrive/Documents (OneDrive)/PHD/RESEARCH/BIOCRUSTS/biocrust_taxtab_10May2023.rds")

