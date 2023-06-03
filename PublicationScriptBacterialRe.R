#Reseqeucning Bacteria
#Preamble for R 
rm(list = ls())
library(ShortRead)
library(dada2)
library(RColorBrewer)
library(phyloseq)
library(ggplot2)
library(vegan)
library(kmer)
library(DESeq2)
#Set working directory
setwd("<directory>")
#Set working directory for trimmed reads.
rawDir <- "./Raw"

#check files list and make lists of the forward sequences ONLY
# get sample names
list.files(rawDir)
fnFs <- sort(list.files(rawDir, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(rawDir, pattern="_R2_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "-110220"), `[`, 1)

#Filter and trim
#create filepaths for results
filtFs <- file.path("./Filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("./Filtered", paste0(sample.names, "_R_filt.fastq.gz"))
#trim and filter
#Trimming parameters
#Left trim selected to remove all the taxa primer sequence (increased from 20 to 25 for forward becasue fastqc shows weirdness)
#truncLen selected to truncate to the point that the mean average read phred drops below 30 (roughly)
#all reads will be of length truncLen - trimLeft
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(17, 21),
                     maxN=0, maxEE = 2, truncQ=10, rm.phix=TRUE, truncLen = c(280, 210),
                     compress=TRUE, multithread=10, verbose = TRUE)
perc.reads.retained <- (out[,2] / out[,1]) * 100
out <- cbind(out, perc.reads.retained)

#learn error rates for algorithm and plot
errF <- learnErrors(filtFs, multithread=10, verbose = 1)
errR <- learnErrors(filtRs, multithread=10, verbose = 1)
errFplot <- plotErrors(errF)
errRplot <- plotErrors(errR)

#dereplicate reads
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#apply core sample inference algorithm 
#lots of diagnostics within this "dada" class object
dadaFs <- dada(derepFs, err=errF, multithread=10)
dadaRs <- dada(derepRs, err=errR, multithread=10)

#and merge
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#make amplicon sequence variant table
seqtab <- makeSequenceTable(mergers)

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=10, verbose=TRUE)

#####################################################################################################################
#assign taxonomy
#species level assignemtn can be made for 16s data (not for fungi) but not sure I actually  need that
#"DECIPHER" libirary  and method is  an alternative to this but not done here.
taxa <- assignTaxonomy(seqtab.nochim,
                       "../References/silva_nr_v132_train_set.fa",
                       multithread = 2, verbose = TRUE)


#Rarefy bacterial reads
######################################
#Remove any sequences that dont have the same start position as the majority of reads by matching forst 4 bases
#Repeat for reverse sequences if applicable
#Difficult as there's a fair amount of variability although one "tag" makes up 75% of all reads
#To be cautious will remove all "tags" that make up less than 1% of all sequences
sequences <- colnames(seqtab.nochim)
tags <- substr(sequences, 1, 4)
unqtags <- unique(tags)
tagnums <- c()
tagabun <- c()
for(i in 1:length(unqtags)){
  this <- sum(tags == unqtags[i])
  nthis <- sum(seqtab.nochim[, which(tags == unqtags[i])])
  tagnums <- c(tagnums, this)
  tagabun <- c(tagabun, nthis)
}
names(tagnums) <- unqtags
names(tagabun) <- unqtags
tagnums
tagabun
perc <- (tagabun / sum(tagabun)) * 100
perc[perc >= 1]
#Remove anything that doesn't start with TGGG, TAGG, TCGG, TAAG, TGAG, CAAG, TCGA, TTAG (to match the previous experiment)
#Mostly lines up woth the more than 1% abundance here too
to.keep <- sequences[which(substr(sequences, 1, 4) == "TGGG" |
                             substr(sequences, 1, 4) == "TAGG" |
                             substr(sequences, 1, 4) == "TCGG" | 
                             substr(sequences, 1, 4) == "TAAG" |
                             substr(sequences, 1, 4) == "TGAG" |
                             substr(sequences, 1, 4) == "CAAG" |
                             substr(sequences, 1, 4) == "TCGA" |
                             substr(sequences, 1, 4) == "TTAG") ]
length(sequences)
length(to.keep)
seqtab.nochim <- seqtab.nochim[, which(colnames(seqtab.nochim) %in% to.keep)]
nchar(colnames(seqtab.nochim))

#Repeat fro reverses
sequences <- colnames(seqtab.nochim)
tags <- substr(sequences, nchar(sequences) - 3, nchar(sequences))
unqtags <- unique(tags)
tagnums <- c()
tagabun <- c()
for(i in 1:length(unqtags)){
  this <- sum(tags == unqtags[i])
  nthis <- sum(seqtab.nochim[, which(tags == unqtags[i])])
  tagnums <- c(tagnums, this)
  tagabun <- c(tagabun, nthis)
}
names(tagnums) <- unqtags
names(tagabun) <- unqtags
tagnums
tagabun
perc <- (tagabun / sum(tagabun)) * 100
perc[perc >= 1]
#Remove tags making up more than 1% of abundance of all reads
to.keep <- sequences[which(substr(sequences,nchar(sequences) - 3, nchar(sequences)) == "AACA" |
                             substr(sequences,nchar(sequences) - 3, nchar(sequences)) == "AACG" |
                             substr(sequences,nchar(sequences) - 3, nchar(sequences)) == "AACC") ]
length(sequences)
length(to.keep)
seqtab.nochim <- seqtab.nochim[, which(colnames(seqtab.nochim) %in% to.keep)]
nchar(colnames(seqtab.nochim))

#remove anything defined as ensifer and try to combine any bio reps that would cause loss at rarefication
taxa <- taxa[which(taxa[, "Genus"] != "Ensifer"), ]
seqtab.nochim <- seqtab.nochim[, which(colnames(seqtab.nochim) %in% rownames(taxa))]

seqtab2 <- seqtab.nochim
seqtab2["53e-1-B", ] <- seqtab2["53e-1-B", ] + seqtab2["53e-3-B", ]
seqtab2 <- seqtab2[which(rownames(seqtab2) != "53e-3-B"), ]
seqtab2 <- seqtab2[which(rownames(seqtab2) != "57e-3-B"), ]
seqtab.nochim <- seqtab2
sample.info <- sample.info[which(rownames(sample.info) %in% rownames(seqtab.nochim)), ]

###Bacteria
############################################################################################################
##############################################################################################
##Alpha diversities
#make phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(sample.info))

#calculate various alpha diversity measures included
alpha_diversities <- estimate_richness(ps, split = TRUE)

#Shannon diversity
fit <- aov(Shannon ~ Treatment + Soil.type, data = alpha_diversities)
summary(fit)
fiti <- aov(Shannon ~ Treatment * Soil.type, data = alpha_diversities)
summary(fiti)
TukeyHSD(fiti)

##############################################################################################
#beta diversities and permanovas
test <- adonis(phyloseq::distance(ps, method = "bray") ~ Soil.type + Treatment + Bio.rep, 
               data = data.frame(sample_data(ps)), method  = "bray", permutations = 100)

#transform sample counts
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
#ordinate
ord.nmds.bray <- ordinate(ps.prop, method="MDS", distance="bray")
#percentage explained variation for each axis is in
ord.nmds.bray$values[1, "Relative_eig"]
ord.nmds.bray$values[2, "Relative_eig"]
ord.nmds.bray$values[3, "Relative_eig"]
ord.nmds.bray$values[4, "Relative_eig"]


#constrained ordination
ordered.info <- sample.info
ordered.info$number <- substr(sample.info$Sample.name, 1, 2)
ordered.info <- ordered.info[order(ordered.info$number), ]
ordered.info

dbRDA <- capscale(seqtab.nochim ~ Treatment, data = ordered.info, distance = "bray")
CAPS <- summary(dbRDA)$sites
CAPS <- cbind(CAPS, ordered.info)
CAPS

CAPS$Soil.type <- factor(CAPS$Soil.type, levels = c("WG", "WS", "PD"), labels = c("WG", "WS", "SP"))
CAPS$Treatment <- factor(CAPS$Treatment, levels = c("mock", "1021", "419", "1022"), labels = c("Mock", "Sm1021", "WSM419", "WSM1022"))

summary(dbRDA)$concont
summary(dbRDA)$cont
summary(dbRDA)$tot.chi
names(summary(dbRDA))



########################################################################################################################
########################################################################################################################
###Bacteria
############Filtering for stacked bar charts
total.abundance <- colSums(seqtab.nochim)
rel.abundance <- (total.abundance / sum(total.abundance)) * 100
seqtab.small <- seqtab.nochim[, which(rel.abundance > 0.01)]

####################################################################################
#Here create an ASV database for comparisons and standardisation
ASV2 <- colnames(seqtab.nochim)
Names <- paste0("ASV", 1:length(ASV2), ".2")
ASV2 <- data.frame(ASV2, Names)
save(ASV2, file = "./ASV2.R")

#subset samples for each soil type and make ps's again fro barplots
sample.infoWG <- sample.info[which(sample.info$Soil.type == "WG"), ]
sample.infoWS <- sample.info[which(sample.info$Soil.type == "WS"), ]
sample.infoPD <- sample.info[which(sample.info$Soil.type == "PD"), ]
seqtab.nochimWG <- seqtab.nochim[which(rownames(seqtab.nochim) %in% rownames(sample.infoWG)), ]
seqtab.nochimWS <- seqtab.nochim[which(rownames(seqtab.nochim) %in% rownames(sample.infoWS)), ]
seqtab.nochimPD <- seqtab.nochim[which(rownames(seqtab.nochim) %in% rownames(sample.infoPD)), ]
#make phyloseq objects
psWG <- phyloseq(otu_table(seqtab.nochimWG, taxa_are_rows=FALSE), 
                 sample_data(sample.infoWG), 
                 tax_table(taxa))
psWS <- phyloseq(otu_table(seqtab.nochimWS, taxa_are_rows=FALSE), 
                 sample_data(sample.infoWS), 
                 tax_table(taxa))
psPD <- phyloseq(otu_table(seqtab.nochimPD, taxa_are_rows=FALSE), 
                 sample_data(sample.infoPD), 
                 tax_table(taxa))

to.test <- list(psWG, psWS, psPD)
cons <- list(c("Treatment", "1022", "mock"), 
             c("Treatment", "1021", "mock"),
             c("Treatment", "419", "mock"),
             c("Treatment", "1022", "1021"), 
             c("Treatment", "1022", "419"), 
             c("Treatment", "1021", "419"))


############ Deseq fucntion
make.result.df <- function(ps.object, con){
  ds1 <- phyloseq_to_deseq2(ps.object, ~ Treatment)
  dds1 <- estimateSizeFactors(ds1, type = "poscount")
  dds1 <- estimateDispersions(dds1, fitType = "local")
  dds1 <- nbinomWaldTest(dds1)
  res1 <- DESeq2::results(dds1, contrast = con)
  df1 <- as.data.frame(res1)
  df1$con <- paste(con, collapse = "")
  df1$taxon <- rownames(df1)
  rownames(df1) <- paste("ArbOTU", 1:length(rownames(df1)), sep = "")
  df <- df1[which(df1$padj < 0.05), ]
  if(nrow(df) > 0){
    df$Phylum <- "phylum"
    df$Class <- "class"
    df$Order <- "order"
    df$Family <- "family"
    df$Genus <- "genus"
    for(x in 1:length(rownames(df))){
      this.data <- taxa[which(rownames(taxa) == df$taxon[x]), ]
      df$Phylum[x] <- this.data[[2]]
      df$Class[x] <- this.data[[3]]
      df$Order[x] <- this.data[[4]]
      df$Family[x] <- this.data[[5]]
      df$Genus[x] <- this.data[[6]]
      print(x)
    }
    df$Soil <- sample_data(ps.object)$Soil.type[1]
  }
  return(df)
}
#loops for every ps.object in to.test list
final.df <- c()
for(q in 1:length(to.test)){
  for(i in 1:length(cons)){
    this.result <- make.result.df(to.test[[q]], cons[[i]])
    if(nrow(this.result) > 0){
      final.df <- rbind(final.df, this.result)
    }
    print(q)
  }
}

#rename with standardise ASV names
final.df$X <- as.character(final.df$X)
for(i in 1:length(final.df$X)){
  final.df$X[i] <- as.character(ASV2$Names[which(as.character(ASV2$ASV2) == as.character(final.df$taxon[i]))])
}

################### CCA analysis endosphere
soil.data <- read.csv("<ReadInNutrientData>")

soil.data <- soil.data[is.na(soil.data$SeqName) == FALSE, ]
rownames(soil.data) <- soil.data$SeqName

sample.info <- sample.info[which(rownames(sample.info) %in% rownames(soil.data)), ]
soil.data <- soil.data[which(rownames(soil.data) %in% rownames(sample.info)), ]
seqtab.nochim <- seqtab.nochim[which(rownames(seqtab.nochim) %in% rownames(soil.data)), ]
sample.info$Treatment <- as.character(sample.info$Treatment)
sample.info$Soil.type <- as.character(sample.info$Soil.type)
#EVERYTHING MUST BE IN THE SAME ORDER
soil.data <- soil.data[order(rownames(soil.data)), ]
seqtab.nochim <- seqtab.nochim[order(rownames(seqtab.nochim)), ]
sample.info <- sample.info[order(rownames(sample.info)), ]

###remove the annotations
soil.data <- soil.data[, -c(1:6)]

#calculate z scores to normalise and combine the clusters
soil.data <- scale(soil.data, center = TRUE, scale = TRUE)
soil.data <- as.data.frame(soil.data)
soil.data$Cluster6P <- soil.data$Ag
soil.data$Cluster4P <- (soil.data$Pb + soil.data$Ni + soil.data$Cd + soil.data$Zn + soil.data$Be) / 5
soil.data$Cluster3P <- (soil.data$S + soil.data$P + soil.data$Co + soil.data$Ba + 
                         soil.data$Sr + soil.data$Se + soil.data$K + soil.data$Nitrogen
                       + soil.data$Mn + soil.data$Carbon + soil.data$Mg + soil.data$Ti_sq
                       + soil.data$Cu + soil.data$Ca + soil.data$B) / 15
soil.data$Cluster5P <- (soil.data$Tl + soil.data$Na + soil.data$Rb + soil.data$Cs) / 4
soil.data$Cluster1P <- (soil.data$U + soil.data$Li + soil.data$As + soil.data$Mo) / 4
soil.data$Cluster2P <- (soil.data$Fe + soil.data$Al + soil.data$Cr + soil.data$V) / 4
soil.data <- soil.data[, -(1:33)]
#Fit model with all clusters
cca1 <- cca(seqtab.nochim ~ ., data = soil.data)







