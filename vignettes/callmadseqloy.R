#Load MADseqLOY and processBam functions
source("/PROJECTES/SHARED/EGA_ncarreras/MADseqLOY.R")
source("/PROJECTES/SHARED/EGA_ncarreras/processBam.R")
source("/PROJECTES/SHARED/EGA_ncarreras/checkSex.R")

#Set the path were the files are stored
setwd("/PROJECTES/SHARED/EGA_ncarreras/EGAD00001005339/")
females <- read.table(file="File_names_female.txt")
males <- read.table(file="File_names_male.txt")

files <- c(females$V1,males$V1)
#Select the male samples that has overcome the quality control
metadata <- read.table(file="metadata_qc.txt", sep="\t")
males_qc <- metadata[which(metadata$Gender=="M" & metadata$Quality_control),]

## Llamo la función madseqloy para todos los archivos
sex <- checkSex("./", mc.cores=10)
ans <- madseqloy(metadata$path, mc.cores=10)

save(ans, file = "ans_madseqloy.Rdata")
