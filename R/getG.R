#TODO move some of the description to details? how much of this is from
# QTLseqr? Could we possibly finesse this so that the option we want is an
# appropriate pull request to QTLseqr?

#' @title get G prime stat
#' @description The function is used by \code{\link{runGprimeAnalysis}} to
#' calculate the G statisic G is defined by the equation:
#' \deqn{G = 2*\sum_{i=1}^{4}n_{i}*ln\frac{obs(n_i)}{exp(n_i)}}{G = 2 * \sum n_i * ln(obs(n_i)/exp(n_i))}
#' Where for each SNP, \eqn{n_i} from i = 1 to 4 corresponds to the reference
#' and alternate allele depths for each bulk, as described in the following
#' table: \tabular{rcc}{ Allele \tab High Bulk \tab Low Bulk \cr Reference \tab
#' \eqn{n_1} \tab \eqn{n_2} \cr Alternate \tab \eqn{n_3} \tab \eqn{n_4} \cr}
#' ...and \eqn{obs(n_i)} are the observed allele depths as described in the data
#' frame. Method 1 calculates the G statistic using expected values assuming
#' read depth is equal for all alleles in both bulks: \deqn{exp(n_1) = ((n_1 +
#' n_2)*(n_1 + n_3))/(n_1 + n_2 + n_3 + n_4)} \deqn{exp(n_2) = ((n_2 + n_1)*(n_2
#' + n_4))/(n_1 + n_2 + n_3 + n_4)} etc...
#'
#' @param LowRef A vector of the reference allele depth in the low bulk
#' @param HighRef A vector of the reference allele depth in the high bulk
#' @param LowAlt A vector of the alternate allele depth in the low bulk
#' @param HighAlt A vector of the alternate allele depth in the high bulk
#'
#' @return A vector of G statistic values with the same length as
#'
#' @details DETAILS
#' @seealso \href{https://doi.org/10.1371/journal.pcbi.1002255}{The Statistics
#'   of Bulk Segregant Analysis Using Next Generation Sequencing}
#'   \code{\link{tricubeStat}} for G prime calculation
#'
#' @export
getG <- function(LowRef, HighRef, LowAlt, HighAlt)
{
  exp <- c(
    (LowRef + HighRef) * (LowRef + LowAlt) / (LowRef + HighRef + LowAlt + HighAlt),
    (LowRef + HighRef) * (HighRef + HighAlt) / (LowRef + HighRef + LowAlt + HighAlt),
    (LowRef + LowAlt) * (LowAlt + HighAlt) / (LowRef + HighRef + LowAlt + HighAlt),
    (LowAlt + HighAlt) * (HighRef + HighAlt) / (LowRef + HighRef + LowAlt + HighAlt)
  )
  obs <- c(LowRef, HighRef, LowAlt, HighAlt)

  G <-
    2 * (rowSums(obs * log(
      matrix(obs, ncol = 4) / matrix(exp, ncol = 4)
    )))
  return(G)
}
