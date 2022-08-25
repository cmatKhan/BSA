# TODO complete documetnation of params, add return description

#' @title Analyze stuff
#'
#' @description Now that we have a list with all comparisons for replicates
#'   (object: "replicates"), pools (object: "pools") and all
#'   samples together (allPoolsInOneComparison) we can begin filtering the
#'   samples. First idea is to filter by the following params
#'
#' @param SNPcomparison: lower tenth percentile >
#'   quantile(sum of depth in both bulks, na.rm = T, probs = 0.1))
#' @param maxDepthPercentile depth: higher fifth percentile >
#'   quantile(sum of depth in both bulks, na.rm = T, probs = 0.95)), Default: 0.1
#' @param minDepthPercentile sample depth - Half the min depth, Default: 0.9
#' @param windowSize PARAM_DESCRIPTION, Default: 25000
#' @param bulkSize PARAM_DESCRIPTION, Default: 20
#' @param outlierFilt see \code{\link{QTLseqR::runGprimeAnalysis}}, Default: 'deltaSNP'
#' @param filter_chr_list PARAM_DESCRIPTION, Default: NULL
#'
#' @return OUTPUT_DESCRIPTION
#'
#' @details DETAILS
#'
#' @export
#'
#' @importFrom dplyr filter
#' @importFrom stats quantile
#' @importFrom QTLseqr filterSNPs
analyzer = function(SNPcomparison,
                     minDepthPercentile = 0.1,
                     maxDepthPercentile = 0.9,
                     windowSize = 2.5e4,
                     bulkSize = 20,
                     outlierFilt = "deltaSNP",
                     filter_chr_list = NULL) {

  # filter the SNP comparison frame -- remove na, and filter out
  # chromosomes passed filter_chr_list
  fltr_SNPcomparison = SNPcomparison %>%
    # TODO consider replacing this with complete.cases(.)?
    dplyr::filter(!is.nan(deltaSNP)) %>%
    dplyr::filter(!is.na(SNPindex.LOW) | !is.na(SNPindex.HIGH)) %>%
    # for instance, if you wish to remove chrM, then pass c("chrM") to
    # filter_chr_list
    filter(if(!is.null(filter_chr_list)) !CHROM %in% filter_chr_list else TRUE)

  # calculate percentiles, passed in as minDepthPercentile and
  # maxDepthPercentile on  of DP.LOW and DP.HIGH
  quants <- stats::quantile(fltr_SNPcomparison$DP.LOW + fltr_SNPcomparison$DP.HIGH,
                     na.rm = TRUE,
                     probs = c(minDepthPercentile,
                               maxDepthPercentile))

  df_filtered <- QTLseqr::filterSNPs(
    SNPset = fltr_SNPcomparison,
    minTotalDepth = quants[1],
    maxTotalDepth = quants[2],
    # minSampleDepth = quants[1]/2,  It seems that the libraries from lung and
    # YPD have less depth than the inoculum, so this parameter was messing up
    # the filtering.
    verbose = TRUE
  )

  # note -- this function is slightly modified from QTLseqR
  df_filtered <- runQTLseqAnalysis_local(df_filtered,
    windowSize = windowSize,
    popStruc = "RIL",
    bulkSize = bulkSize,
    replications = 10000, intervals = c(95, 99)
  )

  # note -- this function is slightly modified from QTLseqR
  df_filtered <- runGprimeAnalysis_local(df_filtered,
                                   windowSize = windowSize,
                                   outlierFilter = outlierFilt,
                                   filterThreshold = 0.1)

  return(df_filtered)
}

