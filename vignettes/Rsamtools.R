library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
library(QDNAseq)

#Choose two files to compare the read counts
setwd("/PROJECTES/SHARED/EGA_ncarreras/EGAD00001005339/")
filename1 <- "EGAF00002727133/PGDX18700P_WGS.sorted_processed.bam" #Male
filename2 <- "EGAF00002727136/PGDX18702P_WGS.sorted_processed.bam" #Male
filename3 <- "EGAF00002727394/PGDX16571P_WGS.sorted_processed.bam" #Female
filename4 <- "EGAF00002727395/PGDX16579P_WGS_X1.sorted_processed.bam" #Female

files <- c(filename1,filename2,filename3,filename4)

#Select the ranges for ChrY
subsetA <- GenomicRanges::GRanges("chrY",IRanges(6611498, 24510581))

#Select the ranges for the autosomes
chr <- paste0("chr",1:22)
bf <- BamFile(filename1)
sl <- seqlengths(bf)
subsetB <- GenomicRanges::GRanges(seqnames = chr, ranges = IRanges::IRanges(1, sl[chr]))


##############
## countBam ##
##############

#Count reads per region
subA <- countBam(filename1, param=ScanBamParam(which = subsetA))
subB <- countBam(filename1, param=ScanBamParam(which = subsetB))

#Calculate coverages
covA <- subA$records/subA$width
covB <- sum(subB$records)/sum(subB$width)

#Ratio
covA/covB


#######################
## summarizeOverlaps ##
#######################

soA <- summarizeOverlaps(features=subsetA,
                         reads=bf,
                         ignore.strand=TRUE)

assay(soA)

soB <- summarizeOverlaps(features=subsetB,
                         reads=bf,
                         ignore.strand=TRUE)

assay(soB)


###############
## MADseqLOY ##
###############


res <- processBam(file=filename1, bins=bins, n.bases.min=30,
                  mappability.min=50, blacklist.max=0, correction=c("gc","mappability"), 
                  span=0.75, family="gaussian", normalization="median", target=subsetA, 
                  ref=subsetB, trim.target=0, trim.ref=0.25)
#0.526/1.005

#No filtering, no correction, no normalization
res <- processBam(file=filename1, bins=bins, n.bases.min=0,
                  mappability.min=0, blacklist.max=100, correction="None", 
                  span=0.75, family="gaussian", normalization="None", target=subsetA, 
                  ref=subsetB, trim.target=0, trim.ref=0.25)

#1245.856/2566.782

res <- processBam(file=filename1, bins=bins, n.bases.min=10,
                  mappability.min=70, blacklist.max=30, correction="gc", 
                  span=0.5, family="symmetric", normalization="mean", target=subsetA, 
                  ref=subsetB, trim.target=0.1, trim.ref=0)

#0.4826685/1.028095

source("/PROJECTES/SHARED/EGA_ncarreras/MADseqLOY.R")
source("/PROJECTES/SHARED/EGA_ncarreras/processBam.R")

res <- madseqloy(fns, mc.cores=2)