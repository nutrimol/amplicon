setwd("~/analysis/amplicon")
library(dada2); packageVersion("dada2")
library(ggplot2); packageVersion("ggplot2")
library(plotly); packageVersion("plotly")
library(GUniFrac); packageVersion("GUniFrac")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(tidyverse); packageVersion("tidyverse")

list.files("./raw/cutadapt")

fnFs <- sort(list.files(path="./raw/cutadapt", pattern="1_001.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path="./raw/cutadapt", pattern="2_001.fq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_L"), `[`, 1)
print("Forward Files list")
fnFs
print("Reverse Files list")
fnRs
print("Sample names")
sample.names

#-------------------------------------Inspect read quality profiles-----------------------------
forwplot<-ggplotly(plotQualityProfile(fnFs[1:length(fnFs)], aggregate = TRUE) +
                     geom_hline(yintercept=c(15,25,35), color=c("red","blue","green"), size=0.5), width =750)
forwplot

revqplot<-ggplotly(plotQualityProfile(fnRs[1:length(fnRs)], aggregate=TRUE) + 
                     geom_hline(yintercept=c(15,25,35), color=c("red","blue","green"), size=0.5), width =750 )
revqplot
#------------------------------------Filter and trim---------------------------------------------
filtFs <- file.path("./raw", "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("./raw", "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(255,221),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, verbose=TRUE, matchIDs = TRUE)
#------------------------------------object out---------------------------------------------
out
print("Total Reads")
sum(out[,2])
#------------------------------------filtered reads-----------------------------------------
plotQualityProfile(filtFs[1:length(fnFs)], aggregate = TRUE) 
plotQualityProfile(filtRs[1:length(fnFs)], aggregate = TRUE) 
#-------------------------------Generating an error model of our data-----------------------
errF <- learnErrors(filtFs,multithread=TRUE, nbases=10000000, MAX_CONSIST = 12, verbose = TRUE)
errR <- learnErrors(filtRs, multithread=TRUE, nbases=10000000, MAX_CONSIST = 12, verbose = TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

drp_F <- derepFastq(filtFs, verbose=TRUE)
drp_R <- derepFastq(filtRs, verbose=TRUE)
#-------------------------------Inferring ASV----------------------------------
ddF <- dada(drp_F, err=errF, multithread=TRUE)
ddR <- dada(drp_R, err=errR, multithread=TRUE)
#----Inferring ASVs using priors to increase sensitivity and detect singletons
sq.etF <- "TGGGGAATATTGCACAATGGGCGAAAGCCTGATGCAGCAACGCCGCGTGAGGGATGAAGGCCTTCGGGTCGTAAACCTCTGTCCTTGGGGAAGAAACAAATGACGGTACCCATGGAGGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCGAGCGTTATCCGGAATTATTGGGCGTAAAGAGTGCGTAGGTGGTTACCTAAGCGCAGGGTCTAAGGCAATGGCTCAACCATTGTTCGCC"
sq.etR <- "CGTTCGCTCCCCATGCTTTCGCGCCTCAGCGTCAGGTTCAGTCCAGAAAGTCGCCTTCGCCACTGGTGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCACTTTCCTCTCCTGCCCTCAAGACAGACAGTTTCAAATGCACGCTCCAGGTTGAGCCCGGAGATTTCACATCTGACTTATCCGTCCGCCTACGCGCCCTTTACGCCCAG"
dd.etF <- dada(drp_F, err=errF, multithread=TRUE, priors = sq.etF, verbose=TRUE)
dd.etR <- dada(drp_R, err=errR, multithread=TRUE, priors = sq.etF, verbose=TRUE)
#---------------------REMOVE BIMERAS DE NOVO------------------------------------
dadaF <- removeBimeraDenovo(makeSequenceTable(ddF), multithread=TRUE, verbose=TRUE)
dadaR <- removeBimeraDenovo(makeSequenceTable(ddR), multithread=TRUE, verbose=TRUE)
dada.etF <- removeBimeraDenovo(makeSequenceTable(dd.etF), multithread=TRUE, verbose=TRUE)
dada.etR <- removeBimeraDenovo(makeSequenceTable(dd.etR), multithread=TRUE, verbose=TRUE)
#-----------------Compare the total set of ASVs inferred by each method#--------
identical(sort(getSequences(dadaF)), sort(getSequences(dada.etF)))   #TRUE <- es que encuentran ambos las mismas ASVs
c(naive=sum(dadaF[,sq.etF]>0), prior=sum(dada.etF[,sq.etF]>0))
c(naive=sum(dadaR[,sq.etR]>0), prior=sum(dada.etR[,sq.etR]>0))
# There are more samples in which the E. timonensis sequence variant was detected when the prior was provided, as expected.
df.emergencia <- data.frame(naive=dadaF[,sq.etF], prior=dada.etF[,sq.etF])
df.emergencia$prior.only <- df.emergencia$naive == 0 & df.emergencia$prior > 0
ggplot(data=df.emergencia, aes(x=naive, y=prior, color=prior.only)) + geom_jitter(alpha=0.5, width=0.2, height=0.2) +
  scale_color_manual(values=c("TRUE"="red", "FALSE"="black")) + guides(color=FALSE) + coord_fixed(ratio=1) +
  xlim(0,20) + ylim(0,20) + xlab("E.timonensis Reads (naive)") + ylab("E. timonensis reads (prior)")
# ------------------------ ASVs inference by psudopooling  y prior----------------------
dd.pseudoF <- dada(drp_F, err=errF, multithread=TRUE, priors = sq.etF, pool="pseudo",verbose=TRUE)
dd.pseudoR <- dada(drp_R, err=errR, multithread=TRUE, priors = sq.etR, pool="pseudo",verbose=TRUE)
dada.pseudoF <- removeBimeraDenovo(makeSequenceTable(dd.pseudoF), multithread=TRUE, verbose=TRUE)
dada.pseudoR <- removeBimeraDenovo(makeSequenceTable(dd.pseudoR), multithread=TRUE, verbose=TRUE)
c(naive=sum(dadaF[,sq.etF]>0), prior=sum(dada.pseudoF[,sq.etF]>0))
c(naive=sum(dadaR[,sq.etR]>0), prior=sum(dada.pseudoR[,sq.etR]>0))
# ------------------------ ASVs inference by psudopooling        ----------------------
dd.pseudoonlyF <- dada(drp_F, err=errF, multithread=TRUE, pool="pseudo",verbose=TRUE)
dd.pseudoonlyR <- dada(drp_R, err=errR, multithread=TRUE, pool="pseudo",verbose=TRUE)
dada.pseudoonlyF <- removeBimeraDenovo(makeSequenceTable(dd.pseudoonlyF), multithread=TRUE, verbose=TRUE)
dada.pseudoonlyR <- removeBimeraDenovo(makeSequenceTable(dd.pseudoonlyR), multithread=TRUE, verbose=TRUE)
c(naive=sum(dadaF[,sq.etF]>0), prior=sum(dada.pseudoonlyF[,sq.etF]>0))
c(naive=sum(dadaR[,sq.etR]>0), prior=sum(dada.pseudoonlyR[,sq.etR]>0))
#---------------------------