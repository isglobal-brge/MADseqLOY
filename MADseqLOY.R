#' Check the ratio between target and reference coverage to detect loss of 
#' chromosome Y events based on sequencing data
#' 
#' @param files A folder path where the .bam files are or a vector of file
#'   paths. This function searches in recursive folders.
#' @param binsize A numeric scalar specifying the width of the bins in units 
#'   of kbp (1000 base pairs). Choices are 1, 5, 10, 15, 30, 50, 100 (default),
#'   500, or 1000.
#' @param genome A character string specifying the genome and genome version to 
#'   be used. Choices are "hg19" (default), "hg18", or "GRCh38".
#' @param target.region The chromosome or region to be compared with the other
#'   regions. By default is the region Y:2694521-59034049 but it can be
#'   changed.
#' @param ref.region If declared, the chromosome or region to be compared with 
#'   the Y region in UCSC style (i.e. "chr21" or "chr21:1000-10000"). By default, 
#'   it includes all the autosomes.
#' @param n.bases.min A numeric scalar between 0 and 100 specifying the minimum
#'   percentage of bases characterized (wihout N). By default is 80.
#' @param mappability.min A numeric scalar between 0 (no reads filtered) and 100 
#'   (maximum reads filtered) specifying the minimum mappability to filter bins. 
#'   By default is 50.
#' @param blacklist.max A numeric scalar between 0 (no overlap) and 100 (maximum
#'   overlap) specifying the maximum overlap between bins and blacklist elements.
#'   By default is 0 (no overlapping allowed).
#' @param correction A character vector specifying which variables to include in 
#'   the correction. Can be c("gc", "mappability") (default), "gc", "mappability",
#'   or "None".
#' @param span The parameter alpha which controls the degree of smoothing in the
#'   read counts correction. By default is 0.75.
#' @param family If correction object different to "None", the loess family. if 
#'   "gaussian" (default) fitting is by least-squares, and if "symmetric" a 
#'   re-descending M estimator is used with Tukey's biweight function.
#' @param normalization A character string specifying the normalization method
#'   for the read counts. Choices are "mean", "median" (default), "mode", or
#'   "None".
#' @param trim.target The fraction (0 to 0.5) of observations to be trimmed
#'   when calculating the coverage of the target.region. By default is 0 (not 
#'   trimmed).
#' @param trim.ref The fraction (0 to 0.5) of observations to be trimmed when 
#'   calculating the coverage of the ref.region. By default is 0.25 (it trimms 
#'   25% observations from both ends).
#' @param mc.cores Number of cores to be used with this function. If there are
#'   more cores than samples, the number of cores will be limited to the number
#'   of samples. By default set to 1.
#' @param qc.sds Threshold to perform quality control of samples using the 
#'   reference region/chromosome. If NULL, the cutoff for considering a good 
#'   quality sample is estimated as 2 times the standard deviation of all samples.
#' @param quiet Should the function not inform about the status of the process.
#'   By default is FALSE.
#' @param ... Other parameters of the functions from QDNAseq package.
#' @return A MADseqLOY object that contains the ratio between target and reference
#'   regions for all the files.
#'   
#' @export
#' @examples
#' \dontrun{
#' madseqloy(files=bamFile)}

madseqloy <- function(files, binsize = 100, genome = "hg19", 
                      target.region, ref.region, n.bases.min = 80, 
                      mappability.min = 50, blacklist.max = 0, 
                      correction=c("gc", "mappability"), span = 0.75, 
                      family = "gaussian", normalization = "median",
                      trim.target = 0, trim.ref = 0.25, mc.cores, qc.sds, 
                      quiet = FALSE, ...){
  
  # Check input-----------------------------------------------------------------
  
  if (missing(files)) {
    stop("A vector of .bam files, a single .bam file or a folder path containing .bam files must be provided\n")
  }
  
  else {
    if (length(files) == 1) {
      if (!file.exists(files)) {
        stop("The given path is neither a file or a folder\n")
      }
      else {
        if (file.info(files)$isdir) {
          allfiles <- grep(pattern = ".bam$", x = list.files(files, 
                                                         recursive = T, full.names = T), perl = T, 
                           value = T)
          if (length(allfiles) == 0) 
            stop("There are no files in the given folder\n")
        }
        else {
          if (length(grep(".bam$", files, perl = T, value = T)) == 
              0) 
            stop("The given file is not a .bam file\n")
          allfiles <- files
        }
      }
    }
    else {
      if (any(!file.exists(files))) {
        if (!quiet) 
          message(paste0("The file ", files[!file.exists(files)], 
                         " does not exist\n"))
        files <- files[file.exists(files)]
      }
      if (!(length(grep(".bam$", files, perl = T, value = T)) == 
            length(files))) {
        if (!quiet) 
          message(paste0("The file ", grep(".bam$", files, 
                                           perl = T, value = T, invert = T), " is not a .bam file\n"))
        files <- grep(".bam$", files, perl = T, value = T)
      }
      allfiles <- files
    }
    if (!quiet) 
      message(paste0("Processing ", length(allfiles), " file(s)...\n"))
  }
  
  # Check bin size -------------------------------------------------------------
  
  if (!binsize %in% c(1,5,10,15,30,50,100,500,1000))
    stop("The binsize should be one of the following ones: 1,5,10,15,30,50,100,500, or 1000.")
  
  # Check hg version and name --------------------------------------------------
  
  if (!genome %in% c("hg18", "hg19", "GRCh38"))
    stop("The human genome release in the hg field should be one of the following ones: 'hg18', 'hg19' or 'GRCh38'.")
  
  # Check target region --------------------------------------------------------
  
  if (missing(target.region)) {
    target.region <- "Y:6611498-24510581"
    message(paste0("Targeted region set to ",  target.region, " by default\n"))
  }
  queryA <- unlist(strsplit(x = target.region, split = "[:, -]", perl = T))
  if (is.na(queryA[2]) | is.na(queryA[3])) {
    subsetA <- GenomicRanges::GRanges(seqnames = gsub("chr", "", queryA[1]), ranges = IRanges::IRanges(start = 1, 
                                                                                                       end = sl[queryA[1]]))
  } else {
    subsetA <- GenomicRanges::GRanges(seqnames = gsub("chr", "", queryA[1]), ranges = IRanges::IRanges(start = as.numeric(queryA[2]), 
                                                                                                       end = as.numeric(queryA[3])))
  }
  
  # Check reference region -----------------------------------------------------
  
  if (missing(ref.region)) {
    message("Using all autosomes as reference region\n")
    subsetB <- GenomicRanges::GRanges(seqnames = 1:22, ranges = IRanges::IRanges(0, 3e+08))
  } else {
    queryB <- unlist(strsplit(x = ref.region, split = "[:, -]", perl = T))
    if (is.na(queryB[2]) | is.na(queryB[3])) {
      subsetB <- GenomicRanges::GRanges(seqnames = gsub("chr", "", queryB[1]), ranges = IRanges::IRanges(start = 1, 
                                                                                                         end = 3e+08))
    } else {
      subsetB <- GenomicRanges::GRanges(seqnames = gsub("chr", "", queryB[1]), ranges = IRanges::IRanges(start = as.numeric(queryB[2]), 
                                                                                                         end = as.numeric(queryB[3])))
    }
  }
  
  # Check correction variables -------------------------------------------------
  
  for (i in correction)  if (!i %in% c("gc","mappability"))
    stop("Correction variables should be gc, mappability, or gc and mappability\n")
  
  # Check Normalization method -------------------------------------------------
  
  if (!normalization %in% c("mean","median","mode","None"))
    stop("The normalization method should be one of the following ones: mean, median, mode, or None.")
  
  # Set mc.cores ---------------------------------------------------------------
  
  if (missing(mc.cores)) 
    mc.cores <- 1
  if (mc.cores > length(allfiles)) {
    mc.cores <- length(allfiles)
    if (!quiet) 
      message(paste0("There are more cores than files to be processed. Parameter 'mc.cores' set to ", 
                     mc.cores, "\n"))
  }
  
  # Obtain reference bins ------------------------------------------------------
  
  bins <- QDNAseq::getBinAnnotations(binSize=binsize, genome=genome)
  
  # Calculate coverage for both regions----------------------------------------

  coverage <- parallel::mclapply(X = allfiles, FUN = processBam, bins = bins, 
                                 n.bases.min = n.bases.min,
                                 mappability.min = mappability.min,
                                 blacklist.max = blacklist.max,
                                 correction = correction, span = span, 
                                 family = family,
                                 normalization = normalization,
                                 target = subsetA, ref = subsetB,
                                 trim.target = trim.target, trim.ref = trim.ref,
                                 mc.cores = mc.cores)
  names(coverage) <- basename(allfiles)
    

  ratio <- unlist(lapply(coverage, "[[", "cov.target")) /
    unlist(lapply(coverage, "[[", "cov.ref"))
  
  # Quality control: standard deviation in reference regions -------------------
  
  sds <- sapply(coverage, "[[", "sd.ref")
  
  if (missing(qc.sds)){
    qc.sds <- 2 * mean(sds)
  }  
  
  ref.qc <- sds > qc.sds
  ratio[ref.qc] <- NA
  
  # Generate output data -------------------------------------------------------
  
  par <- list(files = basename(allfiles),
              hg = genome,
              path = dirname(allfiles),
              target.region = subsetA, 
              ref.region = subsetB,
              trim.ref = trim.ref,
              trim.target = trim.target,
              QCremoved = ref.qc)
  
  Covsummary <- list(Ratio = ratio,
                     Coverage = coverage, 
                     par = par)
  
  class(Covsummary) <- "MADseqLOY"
  return(Covsummary)
}
