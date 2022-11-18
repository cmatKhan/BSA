#' Identify QTL using a smoothed G statistic
#'
#' @note THIS IS NEARLY A VERBATIM COPY OF QTLseqR
#'
#' @description A wrapper for all the functions that perform the full G prime analysis to
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
#' @param SNPset Data frame SNP set containing previously filtered SNPs
#' @param windowSize the window size (in base pairs) bracketing each SNP for which
#'   to calculate the statitics. Magwene et. al recommend a window size of ~25
#'   cM, but also recommend optionally trying several window sizes to test if
#'   peaks are over- or undersmoothed.
#' @param popStruc the population structure. Defaults to "F2" and assumes "RIL" otherwise.
#' @param bulkSize non-negative integer vector. The number of individuals in
#'   each simulated bulk. Can be of length 1, then both bulks are set to the
#'   same size. Assumes the first value in the vector is the simulated high
#'   bulk.
#' @param replications integer. The number of bootstrap replications.
#' @param filter numeric. A minimum SNP-index filter
#' @param intervals confidence intervals -- note this is part of daniel's
#'   additional code to the QTLseqR package
#' @param chrom_col name of the column which stores the chromosome label.
#'   Default is CHROM
#' @param coordinate_col name of the column which stores the snp coordinate.
#'   Default is POS
#' @param delta_snp_col name of the column which stores the change in allele
#'   frequency. Default is deltaSNP
#' @param dp_low_col name of the column which stores the dp_low value. Default
#'   DP.LOW
#' @param dp_high_col name of the column which stores the dp_high value.
#'   Default is DP.HIGH
#'
#' @return The supplied SNP set tibble after G' analysis. Includes five new
#'   columns: \itemize{\item{G - The G statistic for each SNP} \item{Gprime -
#'   The tricube smoothed G statistic based on the supplied window size}
#'   \item{pvalue - the pvalue at each SNP calculatd by non-parametric
#'   estimation} \item{negLog10Pval - the -Log10(pvalue) supplied for quick
#'   plotting} \item{qvalue - the Benajamini-Hochberg adjusted p-value}}
#'
#' @details DETAILS
#'
#' @export
#'
#' @importFrom dplyr group_by mutate left_join
#' @importFrom QTLseqr countSNPs_cpp simulateConfInt
runQTLseqAnalysis_local <- function(SNPset,
                                    windowSize = 1e6,
                                    popStruc = "F2",
                                    bulkSize,
                                    depth = NULL,
                                    replications = 10000,
                                    filter = 0.3,
                                    intervals = c(95, 99),
                                    chrom_col = "CHROM",
                                    coordinate_col = "POS",
                                    delta_snp_col = "deltaSNP",
                                    dp_low_col = "DP.LOW",
                                    dp_high_col = "DP.HIGH") {

  expected_columns = c(chrom_col, coordinate_col, delta_snp_col,
                       dp_low_col, dp_high_col)

  if (length(setdiff(expected_columns, colnames(SNPset))) > 0){
    stop(sprintf("One of %s not in the colnames of input SNPset frame",
                 paste(expected_columns, collapse = ",")))
  }

  message("Counting SNPs in each window...")
  SNPset <- SNPset %>%
    dplyr::group_by(!!rlang::sym(chrom_col)) %>%
    dplyr::mutate(nSNPs =
                    QTLseqr::countSNPs_cpp(
                      POS = !!rlang::sym(coordinate_col),
                      windowSize = windowSize))

  message("Calculating tricube smoothed delta SNP index...")
  SNPset <- SNPset %>%
    dplyr::mutate(tricubeDeltaSNP =
                    tricubeStat_local(POS = !!rlang::sym(coordinate_col),
                                      Stat = !!rlang::sym(delta_snp_col),
                                      windowSize))

  #convert intervals to quantiles
  if (all(intervals >= 1)) {
    message("Returning the following two sided confidence intervals: ",
            paste(intervals, collapse = ", "))
    quantiles <- (100 - intervals) / 200
  } else {
    stop(
      paste0("Convidence intervals ('intervals' paramater) should be ",
             "supplied as two-sided percentiles. i.e. If intervals = '95' will ",
             "return the two sided 95% confidence interval, 2.5% on each side.")
    )
  }

  # calculate min depth per snp between bulks
  SNPset <-
    SNPset %>%
    dplyr::mutate(minDP = pmin(!!rlang::sym(dp_low_col),
                               !!rlang::sym(dp_high_col)))
  SNPset <-
    SNPset %>%
    dplyr::group_by(!!rlang::sym(chrom_col)) %>%
    dplyr::mutate(tricubeDP =
                    floor(tricubeStat_local(!!rlang::sym(coordinate_col),
                                            # TODO figure out where minDP comes
                                            # from
                                            minDP,
                                            windowSize = windowSize)))


  if (is.null(depth)) {
    message(
      "Variable 'depth' not defined, using min and max depth from data: ",
      # TODO figure out where minDP comes from
      min(SNPset$minDP),
      "-",
      # TODO see above
      max(SNPset$minDP)
    )
    # TODO see above
    depth <- min(SNPset$minDP):max(SNPset$minDP)
  }

  #simulate confidence intervals
  CI = simulateConfInt_local(SNPset,
                    popStruc = popStruc,
                    bulkSize = bulkSize,
                    depth = depth,
                    replications = replications,
                    filter = filter,
                    intervals = quantiles
    )


  # match name of column for easier joining of repeat columns
  names(CI)[1] <- "tricubeDP"

  #use join as a quick way to match min depth to matching conf intervals.
  SNPset <-
    dplyr::left_join(x = SNPset,
                     y = CI #, commented out because of above change. need to remove eventually
                     # by = c("tricubeDP" = "depth")
    )

  as.data.frame(SNPset)

}
