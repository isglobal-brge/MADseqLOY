#' Detection of mosaic Loss of chromosome Y in MADseqLOY data
#' 
#' 
#' @param x A MADseqLOY object.
#' @param coef this determines how far the outliers are considered.
#' @param ... Other parameters.
#' @return An object of class 'LOY' that summarizes the LOY events detected in
#'   the analyzed samples
#'   
#' @export
#' @examples
#' \dontrun{
#' getLOY(resMADseqloy)}

getLOY <- function (MADseqLOY, coef=3, ...) {

  if (!inherits(MADseqLOY, "MADseqLOY"))
    stop("MADseqLOY must be an object of class 'MADseqLOY'")

  norm.ratio <- MADseqLOY$Ratio
  
  findOutliers <- function(x, coef=coef) {
    upperq <- stats::quantile(x, 0.75, na.rm=TRUE)
    lowerq <- stats::quantile(x, 0.25, na.rm=TRUE)
    iqr <- stats::IQR(x, na.rm=TRUE)
    # we identify extreme outliers
    extreme.threshold.upper <- (iqr * coef) + upperq
    extreme.threshold.lower <- lowerq - (iqr * coef)
    ans <- x > extreme.threshold.upper | 
      x < extreme.threshold.lower
    ans
  }
  
  getThreshold <- function(x, ...)
  {
    den <- stats::density(x[!is.na(x)], bw="SJ", ...)
    m <- den$x[which(den$y==max(den$y))]
    xx <- x[x>=m]
    x.99 <- stats::quantile(xx, 0.99, na.rm=TRUE)
    ans <- m - x.99
    ans
  }
  
  outl <- findOutliers(norm.ratio, coef=coef)
  
  cl <- rep(NA, length(norm.ratio))
  cl[outl & norm.ratio>0.5] <- "XYY"
  cl[outl & norm.ratio<0.5] <- "LOY"
  cl[!outl] <- "normal"
  cl.f <- factor(cl, levels=c("normal", "LOY", "GOY"))
  
  tt <- getThreshold(norm.ratio)
  fosb <- stats::relevel(cut(norm.ratio, c(-Inf, tt, Inf),
                             labels=c("LOY", "normal")),2)
  ans <- list()
  ans$res <- data.frame(MADseqLOY = cl.f, Fosberg = fosb,
                        continous = norm.ratio)
  ans$data <- MADseqLOY$Ratio
  ans$par <- MADseqLOY$par
  class(ans) <- "LOY"
  
  ans
}
  