#Reseqeucning Fungi
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
#truncLen selected to truncate to similar level to Chrysas run.... merge no problem so nothing gained by having more, retain more reads with shorter truncLen
#all reads will be of length truncLen - trimLeft
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(21, 20),
                     maxN=0, maxEE = 2, truncQ=10, rm.phix=TRUE, truncLen = c(250, 200),
                     compress=TRUE, multithread=10, verbose = TRUE)
perc.reads.retained <- (out[,2] / out[,1]) * 100
out <- cbind(out, perc.reads.retained)
#multithread seems not to make a huge difference
#write out a summary of filtering results
save(out, file = "./Filtered/out.R")

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
taxa <- assignTaxonomy(seqtab.nochim,
                       "../References/UNITE_ITS_v8",
                       multithread = 2, verbose = TRUE)

#####################################################################################################################
#Rarefy fungal reads
sample.info <- read.csv("Sample.info.csv")
rownames(sample.info) <- sample.info[, 1]

#subtract reads in control samples
controls <- seqtab.nochim[37, ]
prop.controls <- controls / sum(controls)
seqtab.nochim <- seqtab.nochim[1:36, ]
for(i in 1:length(rownames(seqtab.nochim))){
  seqtab.nochim[i, ] <- seqtab.nochim[i, ] - (sum(seqtab.nochim[i, ]) * prop.controls)
  print(i)
}
sample.info <- sample.info[which(rownames(sample.info) %in% rownames(seqtab.nochim)), ]
seqtab.nochim[which(seqtab.nochim < 0)] <- 0
seqtab.nochim <- round(seqtab.nochim)

#Need to remove reads associated with a plant.
#Start by identifying the sequences that only classify to the phylum level
library(rBLAST)
potential.plant <-  rownames(taxa)[is.na(taxa[, "Class"])]
potential.plant <- DNAStringSet(potential.plant)
makeblastdb(file = "FILEPATHTOGENOME", dbtype = "nucl")

#Then load library and search
bl <- blast(db = "../Blast/JCVI.Medtr.v4.20130313.fasta")
cl <- predict(bl, potential.plant[1,])
blast.results <- c()
for(i in 1:length(potential.plant)){
  cl <- predict(bl, potential.plant[i, ])
  this.result <- cl[1, ]
  blast.results <- rbind(blast.results, this.result)
  print(i)
}
#create a metric to decide which to dispose of 
#Based on disposing sequences that have at least 90% base identity match over at least 90% of the query length.
my.metric <- (blast.results$Alignment.Length / nchar(potential.plant)) * (blast.results$Perc.Ident / 100)
to.dispose <- which(my.metric >= (0.9*0.9) & blast.results$Perc.Ident >= 90)
#dispose of medicago squences
potential.plant <- potential.plant[to.dispose]
taxa <- taxa[-which(rownames(taxa) %in% potential.plant), ]
seqtab.nochim <- seqtab.nochim[, which(colnames(seqtab.nochim) %in% rownames(taxa))]

#Remove any sequences that dont have the same start position as the majority of reads by matching first 4 bases
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
to.keep <- sequences[which(substr(sequences, 1, 4) == "CGCA") ]
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
to.keep <- sequences[which(substr(sequences,nchar(sequences) - 3, nchar(sequences)) == "TTAA") ]
length(sequences)
length(to.keep)
seqtab.nochim <- seqtab.nochim[, which(colnames(seqtab.nochim) %in% to.keep)]
nchar(colnames(seqtab.nochim))

#rarefy
seqtab.nochim <- rrarefy(seqtab.nochim, sample = 9116)
#remove ASV's with no counts associated now
totals2 <- colSums(seqtab.nochim)
seqtab.nochim <- seqtab.nochim[, which(totals2 != 0)]
#remove singletons (reads that only occur once in one sample)
seqtab.nochim <- seqtab.nochim[, colSums(seqtab.nochim) != 1]

#####################################################################################################################
###Fungi
############Alpha diversities
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

#####################################################################################################################
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

###Fungi
############ Filtering for stacked bar plots
total.abundance <- colSums(seqtab.nochim)
rel.abundance <- (total.abundance / sum(total.abundance)) * 100
seqtab.small <- seqtab.nochim[, which(rel.abundance > 0.01)]

#####################################################################################################################
#Differential abundance analysis 
#Here create an ASV database for comparisons and standardisation
ASV2 <- colnames(seqtab.nochim)
Names <- paste0("ASV", 1:length(ASV2), ".2")
ASV2 <- data.frame(ASV2, Names)

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
  dds1 <- DESeq2::DESeq(ds1)
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

#####################################################################################################################
############################ CCA analysis endosphere

#Read in the nutrient data, soil or plant data
soil.data <- read.csv("<Nutrientdatafile>")
#EVERYTHING MUST BE IN THE SAME ORDER
soil.data <- soil.data[order(rownames(soil.data)), ]
seqtab.nochim <- seqtab.nochim[order(rownames(seqtab.nochim)), ]
sample.info <- sample.info[order(rownames(sample.info)), ]

#calculate z scores to normalise and combine the clusters
soil.data <- scale(soil.data, center = TRUE, scale = TRUE)
soil.data <- as.data.frame(soil.data)
soil.data$Cluster6 <- soil.data$Ag
soil.data$Cluster4 <- (soil.data$Pb + soil.data$Ni + soil.data$Cd + soil.data$Zn + soil.data$Be) / 5
soil.data$Cluster3 <- (soil.data$S + soil.data$P + soil.data$Co + soil.data$Ba + 
                         soil.data$Sr + soil.data$Se + soil.data$K + soil.data$Nitrogen
                       + soil.data$Mn + soil.data$Carbon + soil.data$Mg + soil.data$Ti_sq
                       + soil.data$Cu + soil.data$Ca + soil.data$B) / 15
soil.data$Cluster5 <- (soil.data$Tl + soil.data$Na + soil.data$Rb + soil.data$Cs) / 4
soil.data$Cluster1 <- (soil.data$U + soil.data$Li + soil.data$As + soil.data$Mo) / 4
soil.data$Cluster2 <- (soil.data$Fe + soil.data$Al + soil.data$Cr + soil.data$V) / 4
soil.data <- soil.data[, -(1:33)]
#Fit model with all clusters
cca1 <- cca(seqtab.nochim ~ ., data = soil.data)

#Test the fit of the clusters and axes
cca1.by.term <- anova(cca1, by = "term", perm = 1000)
cca1.by.axis <- anova(cca1, by = "axis", perm = 1000)






