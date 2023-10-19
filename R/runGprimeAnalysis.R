#' @title QTLseqr like G prime analysis
#' @description Identify QTL using a smoothed G statistic
#'
#' @param SNPset Data frame SNP set containing previously filtered SNPs
#' @param windowSize the window size (in base pairs) bracketing each SNP for which
#'   to calculate the statitics. Magwene et. al recommend a window size of ~25
#'   cM, but also recommend optionally trying several window sizes to test if
#'   peaks are over- or undersmoothed. Default: 1e+06
#' @param outlierFilter one of either "deltaSNP" or "Hampel". Method for
#'   filtering outlier (ie QTL) regions for p-value estimation.
#'   Default: 'deltaSNP'
#' @param filterThreshold The absolute delta SNP index to use to filter out putative QTL (default = 0.1)
#' @param ... Other arguments passed to \code{\link[locfit]{locfit}} and
#'   subsequently to \code{\link[locfit]{locfit.raw}}() (or the lfproc). Usefull
#'   in cases where you get "out of vertex space warnings"; Set the maxk higher
#'   than the default 100. See \code{\link[locfit]{locfit.raw}}(). But if you
#'   are getting that warning you should seriously consider increasing your
#'   window size. Default: 0.1
#'
#' @return The supplied SNP set tibble after G' analysis. Includes five new
#'   columns: \itemize{\item{G - The G statistic for each SNP} \item{Gprime -
#'   The tricube smoothed G statistic based on the supplied window size}
#'   \item{pvalue - the pvalue at each SNP calculatd by non-parametric
#'   estimation} \item{negLog10Pval - the -Log10(pvalue) supplied for quick
#'   plotting} \item{qvalue - the Benajamini-Hochberg adjusted p-value}}
#'
#' @details A wrapper for all the functions that perform the full G prime analysis to
#' identify QTL. The following steps are performed:\cr 1) Genome-wide G
#' statistics are calculated by \code{\link{getG}}. \cr G is defined by the
#' equation: \deqn{G = 2*\sum_{i=1}^{4} n_{i}*ln\frac{obs(n_i)}{exp(n_i)}}{G = 2
#' * \sum n_i * ln(obs(n_i)/exp(n_i))} Where for each SNP, \eqn{n_i} from i = 1
#' to 4 corresponds to the reference and alternate allele depths for each bulk,
#' as described in the following table: \tabular{rcc}{ Allele \tab High Bulk
#' \tab Low Bulk \cr Reference \tab \eqn{n_1} \tab \eqn{n_2} \cr Alternate \tab
#' \eqn{n_3} \tab \eqn{n_4} \cr} ...and \eqn{obs(n_i)} are the observed allele
#' depths as described in the data frame. \code{\link{getG}} calculates the G statistic
#' using expected values assuming read depth is equal for all alleles in both
#' bulks: \deqn{exp(n_1) = ((n_1 + n_2)*(n_1 + n_3))/(n_1 + n_2 + n_3 + n_4)}
#' \deqn{exp(n_2) = ((n_2 + n_1)*(n_2 + n_4))/(n_1 + n_2 + n_3 + n_4)}
#' \deqn{exp(n_3) = ((n_3 + n_1)*(n_3 + n_4))/(n_1 + n_2 + n_3 + n_4)}
#' \deqn{exp(n_4) = ((n_4 + n_2)*(n_4 + n_3))/(n_1 + n_2 + n_3 + n_4)}\cr 2) G'
#' - A tricube-smoothed G statistic is predicted by local regression within each
#' chromosome using \code{\link{tricubeStat}}. This works as a weighted average
#' across neighboring SNPs that accounts for Linkage disequilibrium (LD) while
#' minizing noise attributed to SNP calling errors. G values for neighboring
#' SNPs within the window are weighted by physical distance from the focal SNP.
#' \cr \cr 3) P-values are estimated based using the non-parametric method
#' described by Magwene et al. 2011 with the function \code{\link{getPvals}}.
#' Breifly, using the natural log of Gprime a median absolute deviation (MAD) is
#' calculated. The Gprime set is trimmed to exclude outlier regions (i.e. QTL)
#' based on Hampel's rule. An alternate method for filtering out QTL is proposed
#' using absolute delta SNP indeces greater than 0.1 to filter out potential
#' QTL. An estimation of the mode of the trimmed set is calculated using the
#' \code{\link[modeest]{mlv}} function from the package modeest. Finally, the
#' mean and variance of the set are estimated using the median and mode and
#' p-values are estimated from a log normal distribution. \cr \cr 4) Negative
#' Log10- and Benjamini-Hochberg adjusted p-values are calculated using
#' \code{\link[stats]{p.adjust}}
#'
#' @export
#'
#' @importFrom dplyr group_by mutate ungroup
#' @importFrom QTLseqr countSNPs_cpp
#' @importFrom stats p.adjust
runGprimeAnalysis_local <-
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
      dplyr::mutate(tricubeDeltaSNP = tricubeStat_local(POS = POS, Stat = deltaSNP, windowSize, ...))

    message("Calculating G and G' statistics...")
    # note that the warnings from tricuteStat are from the fitting procedure
    # see https://github.com/bmansfeld/QTLseqr/issues/24#issuecomment-521498811
    SNPset %>%
      dplyr::mutate(
        G = getG_local(
          LowRef = AD_REF.LOW,
          HighRef = AD_REF.HIGH,
          LowAlt = AD_ALT.LOW,
          HighAlt = AD_ALT.HIGH)) %>%
      dplyr::mutate(
        Gprime = ifelse(!is.na(G),
                        tricubeStat_local(POS = POS,
                                          Stat = G,
                                          windowSize = windowSize),
                        NA)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        pvalue = ifelse(!is.na(G) & !is.na(Gprime),
                        getPvals_local(
                          Gprime = Gprime,
                          deltaSNP = deltaSNP,
                          outlierFilter = outlierFilter,
                          filterThreshold = filterThreshold),
                        NA)) %>%
      dplyr::mutate(negLog10Pval = -log10(pvalue),
                    qvalue = stats::p.adjust(p = pvalue, method = "BH")) %>%
      as.data.frame()
  }
