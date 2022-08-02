#' @title simulate snp index
#' @description FUNCTION_DESCRIPTION
#'
#' @param depth PARAM_DESCRIPTION
#' @param altFreq1 PARAM_DESCRIPTION
#' @param altFreq2 PARAM_DESCRIPTION
#' @param replicates PARAM_DESCRIPTION, Default: 10000
#' @param filter PARAM_DESCRIPTION, Default: NULL
#'
#' @return OUTPUT_DESCRIPTION
#'
#' @details DETAILS
#'
#' @export
#' @importFrom stats rbinom
simulateSNPindex <-
  function(depth,
           altFreq1,
           altFreq2,
           replicates = 10000,
           filter = NULL) {

    SNPindex_H <- stats::rbinom(replicates, size = depth, altFreq1) / depth
    SNPindex_L <- stats::rbinom(replicates, size = depth, altFreq2) / depth
    deltaSNP <- SNPindex_H - SNPindex_L

    if (!is.null(filter)) {
      deltaSNP <- deltaSNP[SNPindex_H >= filter | SNPindex_L >= filter]
    }
    deltaSNP
  }
