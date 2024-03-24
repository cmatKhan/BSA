#' @export
setClass(
  "BSAExperiment",
  contains = "RangedSummarizedExperiment",
  representation = representation(
    comparisons = "DataFrame"
  )
)

#' @import SummarizedExperiment
BSAExperiment <- function(comparisons = DataFrame(), ...) {

  # Create the SummarizedExperiment object
  se <- SummarizedExperiment(...)

  # Now create the BSAExperiment object, setting comparisons directly
  new("BSAExperiment", se, comparisons = comparisons)
}


S4Vectors::setValidity2("BSAExperiment", function(object) {
  msg <- NULL

  # check that assays has at least DP and AD
  if(is.null(object@assays)){
    warning("`assays` is null", call. = FALSE)
  } else if(!all(c("DP", "AD") %in% names(object@assays))) {
    msg <- c(msg, "The assays must contain at minimum 'DP' and 'AD'")
  }

  # Validate 'comparisons' DataFrame slot
  if (!inherits(object@comparisons, "DataFrame")) {
    msg <- c(msg, "'comparisons' slot must be a DataFrame.")
  }

  if (is.null(msg)) {
    TRUE
  } else msg
})
