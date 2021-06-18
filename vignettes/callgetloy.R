setwd("/PROJECTES/SHARED/EGA_ncarreras/EGAD00001005339/")
load("ans_madseqloy.Rdata")

#Discard samples from quality control and females
metadata <- read.table(file="metadata_qc.txt", sep="\t")
males_qc <- metadata[which(metadata$Gender=="M" & metadata$Quality_control),]

loy <- getLOY(ans)