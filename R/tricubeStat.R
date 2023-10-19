#TODO better title, move some desc into details?

#' @title calc tricube stat
#' @description :ocal regression (wrapper for \code{\link[locfit]{locfit}})
#'   to predict a tricube smoothed version of the statistic supplied for each
#'   SNP. This works as a weighted average across neighboring SNPs that
#'   accounts for Linkage disequilibrium (LD) while minizing noise attributed
#'   to SNP calling errors. Values for neighboring SNPs within the window
#'   are weighted by physical distance from the focal SNP.
#' @param POS A vector of genomic positions for each SNP
#' @param Stat A vector of values for a given statistic for each SNP
#' @param windowSize the window size (in base pairs) bracketing each SNP for which
#'   to calculate the statitics. Magwene et. al recommend a window size of ~25
#'   cM, but also recommend optionally trying several window sizes to test if
#'   peaks are over- or undersmoothed. Default: 2e+06
#' @note example from QTLSeqR: df_filt_4mb$Gprime <- tricubeStat(POS, Stat = GStat, WinSize = 4e6)
#'
#' @return Returns a vector of the weighted statistic caluculted with a
#'   tricube smoothing kernel
#'
#' @details DETAILS
#' @seealso \code{\link{getG}} for G statistic calculation
#' @seealso \code{\link[locfit]{locfit}} for local regression
#'
#' @export
#'
#' @importFrom stats predict
#' @importFrom locfit locfit lp
tricubeStat_local <- function(POS, Stat, windowSize = 2e6)
{
  if (windowSize <= 0)
    stop("A positive smoothing window is required")
  tryCatch(
      stats::predict(
        locfit::locfit(
          Stat ~ locfit::lp(POS, h = windowSize, deg = 0)),
        POS),
    error = function(e){
      browser()
      futile.logger::flog.error(paste0('Error fitting local polynomial model ',
                                'in tricuteStat_local on position: ', POS))
      stop(e)
    })

}
