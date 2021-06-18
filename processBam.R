#' Calculation of the target and reference coverage.
#' 
#' @param file File path of the bam file to be processed.
#' @param bins An AnnotatedDataFrame object with the annotation of all the bins.
#' @param n.bases.min A numeric scalar between 0 and 100 specifying the minimum
#' percentage of bases characterized (wihout N).
#' @param mappability.min A numeric scalar between 0 (no reads filtered) and 100 
#'   (maximum reads filtered) specifying the minimum mappability to filter bins.
#' @param blacklist.max A numeric scalar between 0 (no overlap) and 100 (maximum
#'   overlap) specifying the maximum overlap between bins and blacklist elements.
#' @param correction A character vector specifying which variables to include in 
#'   the correction. Can be c("gc", "mappability"), "gc", "mappability", or 
#'   "None". 
#' @param span If correction object different to "None", the parameter alpha 
#'   which controls the degree of smoothing in the read counts correction.
#' @param family If correction object different to "None", the loess family. if 
#'  "gaussian" fitting is by least-squares, and if "symmetric" a re-descending M 
#'  estimator is used with Tukey's biweight function. 
#' @param normalization A character string specifying the normalization method
#'   for the read counts. Choices are "mean", "median", or "None".
#' @param target A GenomicRanges object with the query region of interest.
#' @param ref A GenomicRanges object with the reference region.
#' @param trim.target The fraction (0 to 0.5) of observations to be trimmed
#'   when calculating the coverage of the target region. 
#' @param trim.ref The fraction (0 to 0.5) of observations to be trimmed
#'   when calculating the coverage of the reference region. 
#' @return A list with the mean and the standard deviation for the reference and
#'   the target region coverages.
#'   
#' @export
#' @examples
#' \dontrun{
#' processBam(file=file.path, bins=getBinAnnotations(100,"hg19"), n.bases.min=30,
#' mappability.min=50, blacklist.max=0, correction=c("gc","mappability"), 
#' span=0.75, family="gaussian", normalization="median", target=GRangestarget, 
#' ref=GRangesref, trim.target=0, trim.ref=0.25)}




processBam <- function(file, bins, n.bases.min, mappability.min, blacklist.max,
                       correction, span, family, normalization, target, ref,
                       trim.target, trim.ref) {
  
  covSummary <- list()
  
  # Count reads per bin---------------------------------------------------------
  
  readCounts <- QDNAseq::binReadCounts(bins, bamfiles=file)
  
  #Transform read counts to GRanges
  binsgr <- GenomicRanges::GRanges(seqnames = bins@data[,"chromosome"], 
                                   ranges = IRanges::IRanges(start = bins@data[,"start"], 
                                                             end = bins@data[,"end"]),
                                   bases = bins@data[,"bases"],
                                   gc = bins@data[,"gc"],
                                   mappability = bins@data[,"mappability"],
                                   blacklist = bins@data[,"blacklist"],
                                   counts = as.numeric(readCounts@assayData$counts))
  
  # Filtering reads-------------------------------------------------------------
  
  message("Filtering reads...\n")
  binsgr_filtered <- binsgr[which(binsgr$bases>=n.bases.min &
                                    binsgr$mappability>=mappability.min &
                                    binsgr$blacklist<=blacklist.max),]
  message(paste0("From a total of ", length(binsgr)," bins, ",length(binsgr)-length(binsgr_filtered)," bins have been filtered.\n"))
  
  # Correcting (loess)----------------------------------------------------------
  
  if (correction[1]=="None"){
    binsgr_filtered$fit <- 1
  }
  else{
    message("Correcting read counts...\n")
    fit <- loess(paste("counts ~", paste(correction, collapse = " * ")), 
                 data=binsgr_filtered, span = span, family = family)
    
    binsgr_filtered$fit <- fit$fitted
  }
  
  binsgr_filtered$counts_norm <- binsgr_filtered$counts/binsgr_filtered$fit
  
  # Normalizing read counts-----------------------------------------------------
  
  if (normalization!="None"){
    message(paste0("Normalizing read counts by ",normalization,"... \n"))
    copynumber <- binsgr_filtered$counts_norm
    if (normalization == "mean") {
      value <- mean(copynumber, na.rm=TRUE)
    } else if (normalization == "median") {
      value <- median(copynumber, na.rm=TRUE)
    } 
    binsgr_filtered$counts_norm <- scale(copynumber, center=FALSE, scale=value)
  }
  
  # Calculate coverage for target region ---------------------------------------
  
  subset.target <- IRanges::subsetByOverlaps(binsgr_filtered,target)
  message("Estimating coverages...\n")
  covSummary$cov.target <- mean(subset.target$counts_norm, na.rm=TRUE, trim=trim.target)
  covSummary$sd.target <- sd(subset.target$counts_norm, na.rm=TRUE)
  
  # Calculate coverage for reference region ------------------------------------
  
  subset.ref <- IRanges::subsetByOverlaps(binsgr_filtered,ref)
  covSummary$cov.ref <- mean(subset.ref$counts_norm, na.rm=TRUE, trim=trim.ref)
  covSummary$sd.ref <- sd(subset.ref$counts_norm, na.rm=TRUE)
  
  return(covSummary)
} 
