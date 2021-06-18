#' Check sex from sequencing data
#' 
#' @param files A folder path where the .bam files are or a vector of file
#'   paths. This function searches in recursive folders.
#' @param mc.cores Number of cores to be used with this function. If there are
#'   more cores than samples, the number of cores will be limited to the number
#'   of samples. By default set to 1.
#' @param quiet Should the function not inform about the status of the process. 
#'   By default is FALSE.
#' @return A checkSex object that contains the trimmed mean coverage in chrX and
#'   chrY in the analyzed samples, and its classification as males and females.
#'   
#' @export
#' @examples
#' \dontrun{
#' checkSex(resMADloy)}

checkSex <- function (files, mc.cores, quiet=FALSE) {
  
  getXY <- function(bamfile) {
    subsetY <- GenomicRanges::GRanges("chrY",IRanges::IRanges(6611498, 24510581))
    subsetX <- GenomicRanges::GRanges("chrX",IRanges::IRanges(2699520, 59034050))
    subsetAuto <- 

    subY <- Rsamtools::countBam(bamfile, 
                                param=Rsamtools::ScanBamParam(which = subsetY))
    subX <- Rsamtools::countBam(bamfile, 
                                param=Rsamtools::ScanBamParam(which = subsetX))

    XYsummary <- list()
    XYsummary$Y <- subY$records/subY$width
    XYsummary$X <- sum(subX$records)/sum(subX$width)

    return(XYsummary)
  }
  
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
  
  # Set mc.cores ---------------------------------------------------------------
  
  if (missing(mc.cores)) 
    mc.cores <- 1
  if (mc.cores > length(allfiles)) {
    mc.cores <- length(allfiles)
    if (!quiet) 
      message(paste0("There are more cores than files to be processed. Parameter 'mc.cores' set to ", 
                     mc.cores, "\n"))
  }
  
  # Get ratio X/Y --------------------------------------------------------------
  
  XY <- parallel::mclapply(allfiles, FUN = getXY, mc.cores = mc.cores)
  data <- do.call(rbind, XY)
  data <- as.data.frame(apply(data, 2, unlist))
  rownames(data) <- basename(allfiles)
  
  # Clustering -----------------------------------------------------------------
  
  pamres <- cluster::pam(data,k=2)
  
  if (pamres$medoids[1]*3 < pamres$medoids[3] & 
     pamres$medoids[2]*3 < pamres$medoids[4]){
    
    class <- factor(rep("FEMALE", length(data$Y)), 
                    levels=c("FEMALE", "MALE"))
    warning("All samples have the same gender status.")
    
  } else if (pamres$medoids[1]*3 > pamres$medoids[3] & 
     pamres$medoids[2]*3 > pamres$medoids[4]){
    
    class <- factor(rep("MALE", length(data$Y)), 
                    levels=c("FEMALE", "MALE"))
    warning("All samples have the same gender status.")
    
  } else if (pamres$medoids[1]*3 < pamres$medoids[3] & 
            pamres$medoids[2]*3 > pamres$medoids[4]){
    
    class <- as.factor(ifelse(pamres$clustering == 1, "FEMALE", "MALE"))
    
  } else if (pamres$medoids[1]*3 > pamres$medoids[3] & 
             pamres$medoids[2]*3 < pamres$medoids[4]){
    
    class <- as.factor(ifelse(pamres$clustering == 2, "FEMALE", "MALE"))
    
  }
  
  par <- list()
  par$path <- files
  par$files <- allfiles
  par$centers <- centers
  
  res <- list(data = data, class = class, par = par)
  class(res) <- "checkSex"
  return(res)
  
}
