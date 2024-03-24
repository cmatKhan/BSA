#' Get the BSAExperiment comparisons frame
#'
#' @rdname BSAExperiment-methods
#' @export
setMethod("comparisons", "BSAExperiment", function(object) {
  object@comparisons
})

#' Set the BSAExperiment comparison frame
#'
#' @rdname BSAExperiment-methods
#'
#' @export
setReplaceMethod("comparisons", "BSAExperiment", function(object, value) {
  # Update object
  object@comparisons <- value
  # validate
  validObject(object)
  # return
  object
})

#' @title Create Sample Comparison Frame
#' @description Create a dataframe which describes comparisons between a given
#'   condition and all other conditions described in a pool dataframe
#' @author chase mateusiak
#'
#' @param pool_df a dataframe which describes the pool
#' @param grouping_var column by which to group and split the
#'   dataframe, Default: 'batch'
#' @param condition_var column which describes the various conditions
#'   of each sample, eg if tissue, then the levels of this column might be
#'   c(lung, brain, YPD, inoculum), Default: 'cond'
#' @param base_comparison_condition condition against which to compare all
#'   other conditions, Default: 'inoculum'
#' @param var1_name the output frame will have two columns, the first storing
#'   the samples which correspond to the `base_comparison_condition`,
#'   the other to the other sample conditions, Default: 'lowBulk'
#' @param var2_name as in var1_name, this will rename the second column
#'   in the output frame, Default: 'highBulk'
#' @param base_cond_in_each_group whether to include the base condition in
#' each group. Default TRUE
#'
#' @return A two column dataframe where the first column is the condition
#'   against which all other sample conditions in that group are compared.
#'   An example structure is:
#'   \tabular{rcc}{ \tab highBulk \tab lowBulk \cr \tab \eqn{P1.1I}
#'   \tab \eqn{P1.1L} \cr \tab
#' \eqn{P1.1I} \tab \eqn{P1.1B} \cr}
#' @details This prepares samples for QTLseqr
#' @examples
#' if(interactive()){
#'    library(dplyr)
#'    # NOTE! "P2.1I" is a singleton
#'     sample_example = c("P1.1I","P1.1L","P1.1Y",
#'                   "P1.2B","P1.2I","P1.2L","P1.2Y","P2.1I")
#'     pool_construction = tibble(sample = sample_example) %>%
#'     mutate(group = str_remove(sample, "\\w$")) %>%
#'    mutate(cond = ifelse(str_detect(sample, "P[[:alnum:]].{1,3}I"),'inoculum', NA)) %>%
#'    mutate(cond = ifelse(str_detect(sample, "P[[:alnum:]].{1,3}Y"),'ypd', cond)) %>%
#'    mutate(cond = ifelse(str_detect(sample, "P[[:alnum:]].{1,3}L"),'lung', cond)) %>%
#'    mutate(cond = ifelse(str_detect(sample, "P[[:alnum:]].{1,3}B"),'brain', cond)) %>%
#'    mutate(bulk = ifelse(cond == "inoculum", 'low', 'high'))
#'    sample_comparison_frame(pool_construction)
#'  }
#' @seealso
#'  \code{\link[dplyr]{reexports}}, \code{\link[dplyr]{group_by}},
#'  \code{\link[dplyr]{group_split}}, \code{\link[dplyr]{pull}},
#'  \code{\link[dplyr]{filter}}, \code{\link[dplyr]{distinct}},
#'  \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{rename}}
#'  \code{\link[rlang]{sym}}
#'  \code{\link[purrr]{map}}
#' @rdname sample_comparison_frame
#' @export
#' @importFrom dplyr setdiff group_by group_split pull filter distinct mutate rename
#' @importFrom rlang sym
#' @importFrom purrr map
#' @importFrom S4Vectors DataFrame
setMethod("createComparisonsFrame",
          "BSAExperiment",
          function(object,
                   grouping_variable = "group",
                   bulk_variable = "bulk",
                   bulk_base_condition = "low",
                   var1_name = 'low_bulk',
                   var2_name = 'high_bulk',
                   base_cond_in_each_group = TRUE) {

  # check colnames of dataframe against expected colnames
  expected_colnames = c(grouping_var, condition_var, "sample")

  if(any(!expected_colnames %in% colnames(pool_df)) & base_cond_in_each_group) {
    stop(paste0("input dataframe colnames must possess a field named sample, ",
                "and at least the grouping_var and conditon_var ",
                sprintf("in any order: %s",
                        paste(expected_colnames, collapse=","))))
  }

  grouping_sym = rlang::sym(grouping_var)
  condition_sym = rlang::sym(condition_var)

  # group the dataframe on grouping_var, and then split into a list of
  # dataframes where each item is one group
  if(base_cond_in_each_group){
    group_split_df = pool_df %>%
      dplyr::group_by(!!grouping_sym) %>%
      dplyr::group_split()

    df_list = purrr::map(group_split_df, function(df) {
      comparison_condition_sample = df %>%
        dplyr::filter(!!condition_sym == base_comparison_condition) %>%
        dplyr::pull(sample)

      expand_grid_df = expand.grid(comparison_condition_sample, df$sample) %>%
        dplyr::distinct(Var2, .keep_all = TRUE) %>%
        dplyr::mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
        dplyr::filter(Var1 != Var2)

      expand_grid_df
    })

    # combine the list of tables into a single table again
    result_df = do.call('rbind', df_list) %>%
      # rename the Var1 and Var2
      dplyr::rename(!!rlang::sym(var1_name) := Var1,
                    !!rlang::sym(var2_name) := Var2)
  } else{
    low_bulk_condition = pool_df %>%
      dplyr::filter(!!condition_sym == base_comparison_condition) %>%
      dplyr::pull(sample)

    if(length(low_bulk_condition) > 1){
      stop("There is more than 1 low bulk condition -- ",
           "setting base_cond_in_each_group failed. Consider splitting up ",
           "the dataframe and running this function on parts, or set ",
           "base_cond_in_each_group to TRUE and allow those conditions which ",
           "do not have base conditions to be dropped.")
    }

    high_bulk_conditions = pool_df %>%
      dplyr::filter(!!condition_sym !=
                      base_comparison_condition) %>%
      dplyr::pull(sample)

    result_df = expand.grid(low_bulk_condition, high_bulk_conditions) %>%
      dplyr::rename(!!rlang::sym(var1_name) := Var1,
                    !!rlang::sym(var2_name) := Var2)
  }

  object@comparisons <- S4Vectors::DataFrame(result_df)

  object
})
