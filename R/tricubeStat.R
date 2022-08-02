#' @title calc tricube stat
#' @description FUNCTION_DESCRIPTION
#'
#' @param POS PARAM_DESCRIPTION
#' @param Stat PARAM_DESCRIPTION
#' @param windowSize PARAM_DESCRIPTION, Default: 2e+06
#'
#' @return OUTPUT_DESCRIPTION
#'
#' @details DETAILS
#'
#' @export
#'
#' @importFrom stats predict
#' @importFrom locfit locfit lp
tricubeStat <- function(POS, Stat, windowSize = 2e6)
{
  if (windowSize <= 0)
    stop("A positive smoothing window is required")
  stats::predict(
    locfit::locfit(
      Stat ~ locfit::lp(POS, h = windowSize, deg = 0)
      ),
    POS)
}
