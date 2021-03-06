% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MADseqLOY.R
\name{madseqloy}
\alias{madseqloy}
\title{Check the ratio between target and reference coverage to detect loss of
chromosome Y events based on sequencing data}
\usage{
madseqloy(
  files,
  binsize = 100,
  genome = "hg19",
  target.region,
  ref.region,
  n.bases.min = 80,
  mappability.min = 50,
  blacklist.max = 0,
  correction = c("gc", "mappability"),
  span = 0.75,
  family = "gaussian",
  normalization = "median",
  trim.target = 0,
  trim.ref = 0.25,
  mc.cores,
  qc.sds,
  quiet = FALSE,
  ...
)
}
\arguments{
\item{files}{A folder path where the .bam files are or a vector of file
paths. This function searches in recursive folders.}

\item{binsize}{A numeric scalar specifying the width of the bins in units
of kbp (1000 base pairs). Choices are 1, 5, 10, 15, 30, 50, 100 (default),
500, or 1000.}

\item{genome}{A character string specifying the genome and genome version to
be used. Choices are "hg19" (default), "hg18", or "GRCh38".}

\item{target.region}{The chromosome or region to be compared with the other
regions. By default is the region Y:2694521-59034049 but it can be
changed.}

\item{ref.region}{If declared, the chromosome or region to be compared with
the Y region in UCSC style (i.e. "chr21" or "chr21:1000-10000"). By default,
it includes all the autosomes.}

\item{n.bases.min}{A numeric scalar between 0 and 100 specifying the minimum
percentage of bases characterized (wihout N). By default is 80.}

\item{mappability.min}{A numeric scalar between 0 (no reads filtered) and 100
(maximum reads filtered) specifying the minimum mappability to filter bins.
By default is 50.}

\item{blacklist.max}{A numeric scalar between 0 (no overlap) and 100 (maximum
overlap) specifying the maximum overlap between bins and blacklist elements.
By default is 0 (no overlapping allowed).}

\item{correction}{A character vector specifying which variables to include in
the correction. Can be c("gc", "mappability") (default), "gc", "mappability",
or "None".}

\item{span}{The parameter alpha which controls the degree of smoothing in the
read counts correction. By default is 0.75.}

\item{family}{If correction object different to "None", the loess family. if
"gaussian" (default) fitting is by least-squares, and if "symmetric" a
re-descending M estimator is used with Tukey's biweight function.}

\item{normalization}{A character string specifying the normalization method
for the read counts. Choices are "mean", "median" (default), "mode", or
"None".}

\item{trim.target}{The fraction (0 to 0.5) of observations to be trimmed
when calculating the coverage of the target.region. By default is 0 (not
trimmed).}

\item{trim.ref}{The fraction (0 to 0.5) of observations to be trimmed when
calculating the coverage of the ref.region. By default is 0.25 (it trimms
25\% observations from both ends).}

\item{mc.cores}{Number of cores to be used with this function. If there are
more cores than samples, the number of cores will be limited to the number
of samples. By default set to 1.}

\item{qc.sds}{Threshold to perform quality control of samples using the
reference region/chromosome. If NULL, the cutoff for considering a good
quality sample is estimated as 2 times the standard deviation of all samples.}

\item{quiet}{Should the function not inform about the status of the process.
By default is FALSE.}

\item{...}{Other parameters of the functions from QDNAseq package.}
}
\value{
A MADseqLOY object that contains the ratio between target and reference
regions for all the files.
}
\description{
Check the ratio between target and reference coverage to detect loss of
chromosome Y events based on sequencing data
}
\examples{
\dontrun{
madseqloy(files=bamFile)}
}
