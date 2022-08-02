#' @title run QTL seq analysis
#' @description FUNCTION_DESCRIPTION
#'
#' @param SNPset PARAM_DESCRIPTION
#' @param windowSize PARAM_DESCRIPTION, Default: 1e+06
#' @param popStruc PARAM_DESCRIPTION, Default: 'F2'
#' @param bulkSize PARAM_DESCRIPTION
#' @param depth PARAM_DESCRIPTION, Default: NULL
#' @param replications PARAM_DESCRIPTION, Default: 10000
#' @param filter PARAM_DESCRIPTION, Default: 0.3
#' @param intervals PARAM_DESCRIPTION, Default: c(95, 99)
#'
#' @return OUTPUT_DESCRIPTION
#'
#' @details DETAILS
#'
#' @export
#'
#' @importFrom dplyr group_by mutate left_join
#' @importFrom QTLseqr countSNPs_cpp simulateConfInt
#' @importFrom timeSeries as.data.frame
runQTLseqAnalysis <- function(SNPset, windowSize = 1e6,
                              popStruc = "F2",
                              bulkSize,
                              depth = NULL,
                              replications = 10000,
                              filter = 0.3,
                              intervals = c(95, 99)
) {

  message("Counting SNPs in each window...")
  SNPset <- SNPset %>%
    dplyr::group_by(CHROM) %>%
    dplyr::mutate(nSNPs = QTLseqr::countSNPs_cpp(POS = POS, windowSize = windowSize))

  message("Calculating tricube smoothed delta SNP index...")
  SNPset <- SNPset %>%
    dplyr::mutate(tricubeDeltaSNP = BSA::tricubeStat(POS = POS, Stat = deltaSNP, windowSize))

  #convert intervals to quantiles
  if (all(intervals >= 1)) {
    message("Returning the following two sided confidence intervals: ", paste(intervals, collapse = ", "))
    quantiles <- (100 - intervals) / 200
  } else {
    stop(
      "Convidence intervals ('intervals' paramater) should be supplied as two-sided percentiles. i.e. If intervals = '95' will return the two sided 95% confidence interval, 2.5% on each side."
    )
  }

  #calculate min depth per snp between bulks
  SNPset <-
    SNPset %>%
    dplyr::mutate(minDP = pmin(DP.LOW, DP.HIGH))

  SNPset <-
    SNPset %>%
    dplyr::group_by(CHROM) %>%
    dplyr::mutate(tricubeDP = floor(BSA::tricubeStat(POS, minDP, windowSize = windowSize)))

  if (is.null(depth)) {
    message(
      "Variable 'depth' not defined, using min and max depth from data: ",
      min(SNPset$minDP),
      "-",
      max(SNPset$minDP)
    )
    depth <- min(SNPset$minDP):max(SNPset$minDP)
  }

  #simulate confidence intervals
  CI <-
    QTLseqr::simulateConfInt(SNPset,
                    popStruc = popStruc,
                    bulkSize = bulkSize,
                    depth = depth,
                    replications = replications,
                    filter = filter,
                    intervals = quantiles
    )


  #match name of column for easier joining of repeat columns
  names(CI)[1] <- "tricubeDP"

  #use join as a quick way to match min depth to matching conf intervals.
  SNPset <-
    dplyr::left_join(x = SNPset,
                     y = CI #, commented out because of above change. need to remove eventually
                     # by = c("tricubeDP" = "depth")
    )

  timeSeries::as.data.frame(SNPset)

}
