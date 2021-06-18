library(ggplot2)

#List files
setwd("/PROJECTES/SHARED/EGA_ncarreras/EGAD00001005339/")
allfiles <- allfiles <- grep(pattern = ".bam$", x = list.files("./", 
                                                               recursive = T, full.names = T), perl = T, 
                             value = T)

#Call getXY
XY <- parallel::mclapply(allfiles, FUN = getXY, mc.cores = 15)
data <- do.call(rbind, XY)
data <- as.data.frame(apply(data, 2, unlist))
rownames(data) <- basename(allfiles)

#Look at the metadata
metadata <- read.table(file="metadata_qc.txt", sep="\t")
rownames(metadata) <- metadata$BAM

metadata_files <- metadata[rownames(data),]
compare <- cbind(metadata_files[,c(7,20)],data)
compare$ratio <- data$ratio

#Clustering
pamres <- cluster::pam(data,k=2)
compare$sex <- pamres$clustering

bam_F <- rownames(compare[which(compare$sex==2&compare$Gender=="F"),])
bam_M <- rownames(compare[which(compare$sex==1&compare$Gender=="M"),])

pamres_males <- cluster::pam(data[bam_M,],k=2)
pamres_females <- cluster::pam(data[bam_F,],k=2)
pamres_mix <- cluster::pam(data[c(bam_F[1:17],bam_M),]$ratio,k=2)

SEX <- as.factor(pamres_mix$clustering)
#Plot results
png("test2.png")
ggplot(data=compare[c(bam_F[1:17],bam_M),],aes(x=1:205,y=ratio,col=SEX))+
  geom_point()
dev.off()

png("test.png")
ggplot(data=compare,aes(x=X,y=Y,col=sex))+
  geom_point()
dev.off()
