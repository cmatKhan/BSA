#TODO move some of description to details? return documentation. better title

#' @title get g prime pvals
#' @description The function is used by \code{\link{runGprimeAnalysis}} to
#' estimate p-values for the weighted G' statistic based on the
#' non-parametric estimation method described in Magwene et al. 2011.
#' Breifly, using the natural log of Gprime a median absolute deviation (MAD)
#' is calculated. The Gprime set is trimmed to exclude outlier regions
#' (i.e. QTL) based on Hampel's rule. An alternate method for filtering out
#' QTL is proposed using absolute delta SNP indeces greater than a set
#' threshold to filter out potential QTL. An estimation of the mode of the
#' trimmed set is calculated using the \code{\link[modeest]{mlv}} function
#' from the package modeest. Finally, the mean and variance of the set are
#' estimated using the median and mode and p-values are estimated from a
#' log normal distribution.
#'
#' @param Gprime a vector of G prime values (tricube weighted G statistics)
#' @param deltaSNP a vector of delta SNP values for use for QTL region filtering
#' @param outlierFilter one of either "deltaSNP" or "Hampel". Method for
#'   filtering outlier (ie QTL) regions for p-value estimation
#' @param filterThreshold The absolute delta SNP index to use to filter out putative QTL
#'
#' @return OUTPUT_DESCRIPTION
#'
#' @details DETAILS
#'
#' @export
#'
#' @importFrom modeest mlv
#' @importFrom stats plnorm
getPvals <-
  function(Gprime,
           deltaSNP = NULL,
           outlierFilter = c("deltaSNP", "Hampel"),
           filterThreshold)
  {

    if (outlierFilter == "deltaSNP") {

      if (abs(filterThreshold) >= 0.5) {
        stop("filterThreshold should be less than 0.5")
      }

      message("Using deltaSNP-index to filter outlier regions with a threshold of ", filterThreshold)
      trimGprime <- Gprime[abs(deltaSNP) < abs(filterThreshold)]
    } else {
      message("Using Hampel's rule to filter outlier regions")
      lnGprime <- log(Gprime)

      medianLogGprime <-median(lnGprime)

      # calculate left median absolute deviation for the trimmed G' prime set
      MAD <-
        median(medianLogGprime - lnGprime[lnGprime <= medianLogGprime])

      # Trim the G prime set to exclude outlier regions (i.e. QTL) using Hampel's rule
      trimGprime <-
        Gprime[lnGprime - median(lnGprime) <= 5.2 * MAD]
    }

    medianTrimGprime <- median(trimGprime)

    # estimate the mode of the trimmed G' prime set using the half-sample method
    message("Estimating the mode of a trimmed G prime set using the 'modeest' package...")
    modeTrimGprime <-
      modeest::mlv(x = trimGprime, bw = 0.5, method = "hsm")[1]

    muE <- log(medianTrimGprime)
    varE <- abs(muE - log(modeTrimGprime))
    #use the log normal distribution to get pvals
    message("Calculating p-values...")
    pval <-
      1 - stats::plnorm(q = Gprime,
                 meanlog = muE,
                 sdlog = sqrt(varE))

    return(pval)
  }
