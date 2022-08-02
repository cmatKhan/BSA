#' @title analyze the g prime
#' @description FUNCTION_DESCRIPTION
#'
#' @param SNPset PARAM_DESCRIPTION
#' @param windowSize PARAM_DESCRIPTION, Default: 1e+06
#' @param outlierFilter PARAM_DESCRIPTION, Default: 'deltaSNP'
#' @param filterThreshold PARAM_DESCRIPTION, Default: 0.1
#' @param ... PARAM_DESCRIPTION
#'
#' @return OUTPUT_DESCRIPTION
#'
#' @details DETAILS
#'
#' @export
#'
#' @importFrom dplyr group_by mutate ungroup
#' @importFrom QTLseqr countSNPs_cpp getPvals
#' @importFrom stats p.adjust
#' @importFrom timeSeries as.data.frame
runGprimeAnalysis <-
  function(SNPset,
           windowSize = 1e6,
           outlierFilter = "deltaSNP",
           filterThreshold = 0.1,
           ...)
  {
    message("Counting SNPs in each window...")
    SNPset <- SNPset %>%
      dplyr::group_by(CHROM) %>%
      dplyr::mutate(nSNPs = QTLseqr::countSNPs_cpp(POS = POS, windowSize = windowSize))

    message("Calculating tricube smoothed delta SNP index...")
    SNPset <- SNPset %>%
      dplyr::mutate(tricubeDeltaSNP = BSA::tricubeStat(POS = POS, Stat = deltaSNP, windowSize, ...))

    message("Calculating G and G' statistics...")
    SNPset <- SNPset %>%
      dplyr::mutate(
        G = BSA::getG(
          LowRef = AD_REF.LOW,
          HighRef = AD_REF.HIGH,
          LowAlt = AD_ALT.LOW,
          HighAlt = AD_ALT.HIGH
        ),
        Gprime = BSA::tricubeStat(
          POS = POS,
          Stat = G,
          windowSize = windowSize,
          ...
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        pvalue = QTLseqr::getPvals(
          Gprime = Gprime,
          deltaSNP = deltaSNP,
          outlierFilter = outlierFilter,
          filterThreshold = filterThreshold
        ),
        negLog10Pval = -log10(pvalue),
        qvalue = stats::p.adjust(p = pvalue, method = "BH")
      )

    return(timeSeries::as.data.frame(SNPset))
  }
