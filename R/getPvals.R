#' @title get g prime pvals
#' @description FUNCTION_DESCRIPTION
#'
#' @param Gprime PARAM_DESCRIPTION
#' @param deltaSNP PARAM_DESCRIPTION, Default: NULL
#' @param outlierFilter PARAM_DESCRIPTION, Default: c("deltaSNP", "Hampel")
#' @param filterThreshold PARAM_DESCRIPTION
#'
#' @return OUTPUT_DESCRIPTION
#'
#' @details DETAILS
#'
#' @export
#'
#' @importFrom timeSeries median
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

      medianLogGprime <- timeSeries::median(lnGprime)

      # calculate left median absolute deviation for the trimmed G' prime set
      MAD <-
        timeSeries::median(medianLogGprime - lnGprime[lnGprime <= medianLogGprime])

      # Trim the G prime set to exclude outlier regions (i.e. QTL) using Hampel's rule
      trimGprime <-
        Gprime[lnGprime - timeSeries::median(lnGprime) <= 5.2 * MAD]
    }

    medianTrimGprime <- timeSeries::median(trimGprime)

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
