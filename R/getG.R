#' @title get G prime stat
#'
#' @description FUNCTION_DESCRIPTION
#'
#' @param LowRef PARAM_DESCRIPTION
#' @param HighRef PARAM_DESCRIPTION
#' @param LowAlt PARAM_DESCRIPTION
#' @param HighAlt PARAM_DESCRIPTION
#'
#' @return OUTPUT_DESCRIPTION
#'
#' @details DETAILS
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
