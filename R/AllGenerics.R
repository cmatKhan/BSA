#' @rdname BSAExperiment-methods
#' @export
setGeneric("comparisons", function(object, comparisons)
  standardGeneric("comparisons")
)

#' @rdname BSAExperiment-methods
#' @export
setGeneric("comparisons<-", function(object, ..., value)
  standardGeneric("comparisons<-")
)

#' @rdname BSAExperiment-methods
#' @export
setGeneric(name = "createComparisonsFrame",
           signature = "object",
           def = function(object,
                          grouping_variable = "batch",
                          bulk_variable = "condition",
                          bulk_base_condition = "inoculum",
                          var1_name = 'low_bulk',
                          var2_name = 'high_bulk',
                          base_cond_in_each_group = TRUE) {
             standardGeneric("createComparisonsFrame")
           })
