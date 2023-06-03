#Master document for creation of fungal sequence tables
#Preamble for R
rm(list = ls())
library(ShortRead)
library(dada2)
library(RColorBrewer)
library(phyloseq)
library(ggplot2)
library(vegan)
library(ape)
library(Biostrings)

#Set working directory
setwd("<directory>")
#Set working directory for trimmed reads.
rawDir <- "./Raw"

#check files list and make lists of the forward sequences ONLY
# get sample names
list.files(rawDir)
fnFs <- sort(list.files(rawDir, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(rawDir, pattern="_R2_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "-271118"), `[`, 1)

#Filter and trim
#create filepaths for results
filtFs <- file.path("./Filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("./Filtered", paste0(sample.names, "_R_filt.fastq.gz"))
#trim and filter
#Trimming parameters
#Left trim selected to remove all the taxa primer sequence (increased from 20 to 25 for forward becasue fastqc shows weirdness)
#truncLen selected to truncate to the point that the mean average read phred drops below 30 (roughly)
#all reads will be of length truncLen - trimLeft
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(21, 20),
                     maxN=0, maxEE = 2, truncQ=10, rm.phix=TRUE, truncLen = c(230, 180),
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
save(errFplot, file = "./ErrorPlots/errFplot.R")
save(errRplot, file = "./ErrorPlots/errRplot.R")
save(errF, file = "./ErrorPlots/errF.R")
save(errR, file = "./ErrorPlots/errR.R")

#dereplicate reads
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
save(derepFs, file = "./Dereplicated/derepFs.R")
save(derepRs, file = "./Dereplicated/derepRs.R")

#apply core sample inference algorithm 
#lots of diagnostics within this "dada" class object
dadaFs <- dada(derepFs, err=errF, multithread=10)
dadaRs <- dada(derepRs, err=errR, multithread=10)
save(dadaFs, file = "./DadaObjects/dadaFs.R")
save(dadaRs, file = "./DadaObjects/dadaRs.R")


#and merge
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
save(mergers, file = "./DadaObjects/mergers.R")

#make amplicon sequence variant table
seqtab <- makeSequenceTable(mergers)
save(seqtab, file = "./SequenceTables/seqtab.R")

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=10, verbose=TRUE)
save(seqtab.nochim, file = "./SequenceTables/seqtab.nochim.R")

#Final track of reads throughout the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out[,1:2], sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
save(track, file = "./track.R")

#####################################################################################################################
#assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim,
                       "../References/UNITE_ITS_v8",
                       multithread = 2, verbose = TRUE)

#inspect assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
#if happy save file out for use
save(taxa, file  = "./taxa.R")


#####################################################################################################################


sample.info <- read.csv("SampleInfo.csv")
rownames(sample.info) <- sample.info$SampleID
sample.info <- sample.info[,3:5]
rownames(seqtab.nochim) <- gsub("-", "_", rownames(seqtab.nochim))


#there are reads in the control samples.
#need to subtract a mean of control reads from all the samples
#Must be carefull here so as to account for each sample's (and control's) total read number
#Take the 4 control samples, rarefy to 500 to encompass all and average reads fro each OTU
#callculate the proportion of reads from each OTU (in controls)
#remove that proportion of reads from each OTU in every sample.
controls <- seqtab.nochim[349:352, ]
controls <- rrarefy(controls, 580)
avg.controls <- colSums(controls) / length(rownames(controls))
prop.controls <- avg.controls / sum(avg.controls)
seqtab.nochim <- seqtab.nochim[1:348, ]
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
#only need to run this onceto  make the blast database for medicago     
#makeblastdb(file = "./AnalysisFungiTax010/Blast/JCVI.Medtr.v4.20130313.fasta", dbtype = "nucl")

#Then load library and search
bl <- blast(db = "<pathtogenome>")
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


#Remove any sequences that dont have the same start position as the majority of reads by matching forst 4 bases
#Repeat for reverse sequences if applicable
sequences <- colnames(seqtab.nochim)
tags <- substr(sequences, 1, 4)
unqtags <- unique(tags)
tagnums <- c()
for(i in 1:length(unqtags)){
  this <- sum(tags == unqtags[i])
  tagnums <- c(tagnums, this)
}
names(tagnums) <- unqtags
tagnums[tagnums != 1]
#Remove anything that doesn't start with CGCA, very clear in this instance

sequences <- colnames(seqtab.nochim)
tags <- substr(sequences, nchar(sequences) - 3, nchar(sequences))
unqtags <- unique(tags)
tagnums <- c()
for(i in 1:length(unqtags)){
  this <- sum(tags == unqtags[i])
  tagnums <- c(tagnums, this)
}
names(tagnums) <- unqtags
tagnums[tagnums != 1]
#Remove anything that doesn't end TTAA

to.keep <- sequences[which(substr(sequences, 1, 4) == "CGCA" & 
                          substr(sequences, nchar(sequences) - 3, nchar(sequences)) == "TTAA")]
length(sequences)
length(to.keep)
seqtab.nochim <- seqtab.nochim[, which(colnames(seqtab.nochim) %in% to.keep)]
nchar(colnames(seqtab.nochim))

#rarefy at 4000 reads
orig.seqtab.nochim <- seqtab.nochim
seqtab.nochim <- rrarefy(orig.seqtab.nochim, sample = 4000)
#remove ASV's with no counts associated now
totals2 <- colSums(seqtab.nochim)
seqtab.nochim <- seqtab.nochim[, which(totals2 != 0)]

#want to merge all the technical replicates into one for each biological replicate and treatment
sample.info$CONCAT <- paste(sample.info$TREATEMENT, 
                               sample.info$FRACTION, 
                               sample.info$SOIL.TYPE, 
                               substr(rownames(sample.info), 1, 1
                               ))
new.samples <- unique(sample.info$CONCAT)
new.table <- matrix(0, length(new.samples), ncol(seqtab.nochim))
rownames(new.table) <- new.samples
for(x in 1:length(new.samples)){
  these.samples <- rownames(sample.info[which(sample.info$CONCAT == new.samples[x]), ])
  if(length(these.samples) > 0){
    if(length(these.samples) > 1){
      sub.table <- seqtab.nochim[which(rownames(seqtab.nochim) %in% these.samples), ]
      rownames(new.table)[x] <- rownames(sub.table)[1]
      new.table[x, ] <- colSums(sub.table) / length(sub.table[, 1])
      print(x)
    }
    else{
      sub.table <- seqtab.nochim[which(rownames(seqtab.nochim) %in% these.samples), ]
      new.table[x, ] <- sub.table
    }
   
  }
}
colnames(new.table) <- colnames(seqtab.nochim)

seqtab.nochim.merged <- new.table
sample.info.merged <- sample.info[which(rownames(sample.info) %in% rownames(seqtab.nochim.merged)), ]
#Have to round too; alpha and betadiversityb calculations dont like non-integers.
seqtab.nochim.merged <- round(seqtab.nochim.merged)
#remove singletons (reads that only occur once in one sample)
seqtab.nochim.merged <- seqtab.nochim.merged[, colSums(seqtab.nochim.merged) > 1]


##############################################################################################################################
#Alpha diversity
#Rep1 problem doesn't exist here for the fungi so remove only the unplanted soils
sample.info.merged <- sample.info.merged[which(sample.info.merged$TREATEMENT != "unplanted"), ]
seqtab.nochim.merged <- seqtab.nochim.merged[which(rownames(seqtab.nochim.merged) %in% rownames(sample.info.merged)), ]

ps <- phyloseq(otu_table(seqtab.nochim.merged, taxa_are_rows=FALSE), 
               sample_data(sample.info.merged))

#calculate various alpha diversity measures included
alpha_diversities <- estimate_richness(ps, split = TRUE)

#add sample information to the data frame
sample.info.merged$TREATEMENT <- as.character(sample.info.merged$TREATEMENT)
sample.info.merged$FRACTION <- as.character(sample.info.merged$FRACTION)
sample.info.merged$SOIL.TYPE <- as.character(sample.info.merged$SOIL.TYPE)
rownames(alpha_diversities) <- gsub("X", "", rownames(alpha_diversities))
alpha_diversities$TREATEMENT <- ""
alpha_diversities$FRACTION <- ""
alpha_diversities$SOIL.TYPE <- ""
for(i in 1:length(rownames(alpha_diversities))){
  alpha_diversities[i, 10:12] <- as.character(sample.info.merged[which(rownames(sample.info.merged) == rownames(alpha_diversities)[i]), 1:3])
}
#add biological rep information
alpha_diversities$Rep <- substr(rownames(alpha_diversities), 1, 1)
alpha_diversities$SOIL.TYPE <- factor(alpha_diversities$SOIL.TYPE, levels=c("WG","WS","PD"), labels=c("WG","WS","PD")) 
alpha_diversities$FRACTION <- factor(alpha_diversities$FRACTION, levels=c("input soil","soil","rhizosphere", "endosphere"), labels=c("Input","Soil","Rhizosphere", "Endosphere")) 

alpha_diversities$combined <- paste(alpha_diversities$TREATEMENT, alpha_diversities$SOIL.TYPE, alpha_diversities$FRACTION)
bartlett.test(Shannon ~ combined, data = alpha_diversities[which(alpha_diversities$FRACTION == "Endosphere"),])
bartlett.test(Shannon ~ combined, data = alpha_diversities[which(alpha_diversities$FRACTION == "Rhizosphere"),])
bartlett.test(Shannon ~ combined, data = alpha_diversities[which(alpha_diversities$FRACTION == "Soil"),])
shapiro.test(alpha_diversities[which(alpha_diversities$FRACTION == "Endosphere"), "Shannon"])
shapiro.test(alpha_diversities[which(alpha_diversities$FRACTION == "Rhizosphere"), "Shannon"])
shapiro.test(alpha_diversities[which(alpha_diversities$FRACTION == "Soil"), "Shannon"])

#Shannon diversity
fit <- aov(Shannon ~ FRACTION + SOIL.TYPE + TREATEMENT, data = alpha_diversities)
fiti <- aov(Shannon ~ FRACTION * SOIL.TYPE * TREATEMENT , data = alpha_diversities)

#Shannon diversity
fite <- aov(Shannon ~ TREATEMENT + SOIL.TYPE, data = alpha_diversities[which(alpha_diversities$FRACTION == "endosphere"), ])
fitr <- aov(Shannon ~ TREATEMENT + SOIL.TYPE, data = alpha_diversities[which(alpha_diversities$FRACTION == "rhizosphere"), ])
fits<- aov(Shannon ~ TREATEMENT + SOIL.TYPE, data = alpha_diversities[which(alpha_diversities$FRACTION == "soil" | alpha_diversities$FRACTION == "input soil"), ])
fitei <- aov(Shannon ~ TREATEMENT * SOIL.TYPE, data = alpha_diversities[which(alpha_diversities$FRACTION == "endosphere"), ])
fitri <- aov(Shannon ~ TREATEMENT * SOIL.TYPE, data = alpha_diversities[which(alpha_diversities$FRACTION == "rhizosphere"), ])
fitsi <- aov(Shannon ~ TREATEMENT * SOIL.TYPE, data = alpha_diversities[which(alpha_diversities$FRACTION == "soil" | alpha_diversities$FRACTION == "input soil"), ])

##############################################################################################

#Rep1 problem doesn't exist here for the fungi so remove only the unplanted soils
sample.info.merged <- sample.info.merged[which(sample.info.merged$TREATEMENT != "unplanted"), ]
seqtab.nochim.merged <- seqtab.nochim.merged[which(rownames(seqtab.nochim.merged) %in% rownames(sample.info.merged)), ]

ps <- phyloseq(otu_table(seqtab.nochim.merged, taxa_are_rows=FALSE), 
               sample_data(sample.info.merged))


#Add biorep infromation where it needs to go
sample_data(ps)$Biorep <- substr(rownames(sample_data(ps)), 1, 1)
#beta diversities and permanovas
test <- adonis(phyloseq::distance(ps, method = "bray") ~ SOIL.TYPE + FRACTION + TREATEMENT + Biorep, 
               data = data.frame(sample_data(ps)), method  = "bray", permutations = 100)
test
capture.output(print("Bray permanova all"), file = "./AnalysisFungiTax011/beta.statistics.txt")
capture.output(test, file = "./AnalysisFungiTax011/beta.statistics.txt", append = TRUE)


#transform sample counts
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
#ordinate
ord.nmds.bray <- ordinate(ps.prop, method="MDS", distance="bray")
#percentage explained variation for each axis is in
ord.nmds.bray$values[1:2, "Relative_eig"]
#make an object for ggplot
vectors <- as.data.frame(ord.nmds.bray$vectors)
vectors <- cbind(vectors, sample_data(ps.prop))

#beta diversities and permanovas for bulk soil only
psbulk <- subset_samples(ps, sample_data(ps)$FRACTION %in% c("soil", "input soil"))
test <- adonis(phyloseq::distance(psbulk, method = "bray") ~ SOIL.TYPE + TREATEMENT + Biorep, 
               data = data.frame(sample_data(psbulk)), method  = "bray", permutations = 100)

#transform sample counts
ps.prop <- transform_sample_counts(psbulk, function(otu) otu/sum(otu))
#ordinate
ord.nmds.bray <- ordinate(ps.prop, method="MDS", distance="bray")
#percentage explained variation for each axis is in
ord.nmds.bray$values[1:2, "Relative_eig"]


#beta diversities and permanovas for rhizosphere only
psrhiz <- subset_samples(ps, sample_data(ps)$FRACTION %in% c("rhizosphere"))
test <- adonis(phyloseq::distance(psrhiz, method = "bray") ~ SOIL.TYPE + TREATEMENT + Biorep, 
               data = data.frame(sample_data(psrhiz)), method  = "bray", permutations = 100)

#transform sample counts
ps.prop <- transform_sample_counts(psrhiz, function(otu) otu/sum(otu))
#ordinate
ord.nmds.bray <- ordinate(ps.prop, method="MDS", distance="bray")
#percentage explained variation for each axis is in
ord.nmds.bray$values[1:2, "Relative_eig"]


#beta diversities and permanovas for endosphere only
psendo <- subset_samples(ps, sample_data(ps)$FRACTION %in% c("endosphere"))
test <- adonis(phyloseq::distance(psendo, method = "bray") ~ SOIL.TYPE + TREATEMENT + Biorep, 
               data = data.frame(sample_data(psendo)), method  = "bray", permutations = 100)

#transform sample counts
ps.prop <- transform_sample_counts(psendo, function(otu) otu/sum(otu))
#ordinate
ord.nmds.bray <- ordinate(ps.prop, method="MDS", distance="bray")
#percentage explained variation for each axis is in
ord.nmds.bray$values[1:2, "Relative_eig"]

############################ CCA analysis all bulk soils

#vettel model building too complex for Mac
soil.data <- read.csv("<NutrientData>")
seqtab.nochim <- seqtab.nochim.merged
sample.info <- sample.info.merged

soil.data <- soil.data[is.na(soil.data$SeqName) == FALSE, ]
rownames(soil.data) <- soil.data$SeqName
sample.info <- sample.info[which(rownames(sample.info) %in% rownames(soil.data)), ]
seqtab.nochim <- seqtab.nochim[which(rownames(seqtab.nochim) %in% rownames(soil.data)), ]
sample.info$TREATEMENT <- as.character(sample.info$TREATEMENT)
sample.info$SOIL.TYPE <- as.character(sample.info$SOIL.TYPE)
#EVERYTHING MUST BE IN THE SAME ORDER
soil.data <- soil.data[order(rownames(soil.data)), ]
seqtab.nochim <- seqtab.nochim[order(rownames(seqtab.nochim)), ]
sample.info <- sample.info[order(rownames(sample.info)), ]

###remove the annotations
soil.data <- soil.data[, -c(1:6)]

#calculate z scores to normalise and combine the clusters
soil.data <- scale(soil.data, center = TRUE, scale = TRUE)
soil.data <- as.data.frame(soil.data)
soil.data$Cluster6S <- (soil.data$Nitrite + soil.data$Pb) / 2
soil.data$Cluster4S <- (soil.data$Rb + soil.data$B + soil.data$Phosphate) / 3
soil.data$Cluster1S <- (soil.data$Ba + soil.data$Al + soil.data$Sr + soil.data$Na + 
                         soil.data$Ca + soil.data$Carbon + soil.data$Ti_sq + soil.data$Cr
                       + soil.data$K + soil.data$Mg + soil.data$P + soil.data$U
                       + soil.data$Tl + soil.data$V + soil.data$Cd + soil.data$Nitrogen + 
                         soil.data$pH + soil.data$Cu + soil.data$As + soil.data$Be) / 20
soil.data$Cluster5S <- (soil.data$Nitrate + soil.data$conductivity + soil.data$Sulphate) / 3
soil.data$Cluster3S <- (soil.data$Ag + soil.data$Ni + soil.data$Zn) / 3
soil.data$Cluster2S <- (soil.data$Fluoride + soil.data$Co + soil.data$Cs + soil.data$Fe + soil.data$Mn + soil.data$Li + soil.data$Mo + soil.data$Se) / 8
soil.data <- soil.data[, -(1:39)]

#Fit model with all clusters
cca1 <- cca(seqtab.nochim ~ ., data = soil.data)


