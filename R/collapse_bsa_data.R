#' Collapse BSA (Bulk Segregant Analysis) Data
#'
#' This function takes a dataframe of BSA data, a list of grouping conditions,
#' and a list of sample columns, and returns a dataframe in which the data has
#' been grouped according to the specified conditions. It performs several
#' operations including grouping, summarizing, uniting, and arranging.
#'
#' @param df_with_meta A dataframe that contains BSA data.
#' @param grouping_conditions A vector of column names that specifies the
#'   grouping conditions.
#' @param sample_columns A vector of column names that specifies the
#'   sample columns.
#'
#' @return A dataframe in which the data has been grouped, summarized,
#'   united, and arranged according to the specified conditions and columns.
#'
#' @importFrom dplyr group_by syms summarize arrange ungroup
#' @importFrom tidyr unite
#'
#' @examples
#' df_with_meta <- data.frame(
#'   A = rep(c("group1", "group2"), each = 5),
#'   CHR = rep("CHR1", 10),
#'   POS = rep(c(1,5,10,20,30), times=2),
#'   REF_Allele = rep('A', 10),
#'   ALT1 = rep('G', 10),
#'   RealDepth = runif(10, 1, 10),
#'   Depth = runif(10, 1, 10),
#'   Reference = runif(10, 1, 10),
#'   Alternative1 = runif(10, 1, 10)
#' )
#' grouping_conditions <- c('CHR','POS','REF_Allele','ALT1')
#' sample_columns <- c("RealDepth", "Depth", "Reference", "Alternative1")
#'
#' collapse_bsa_data(df_with_meta, grouping_conditions, sample_columns)
#'
#' @export
collapse_bsa_data <- function(df_with_meta,
                              grouping_conditions,
                              sample_columns) {
  req_cols = union(grouping_conditions,
               c('CHR', 'POS', 'REF_Allele', 'ALT1',
                 'RealDepth', 'Depth', 'Reference', 'Alternative1'))
  if (length(setdiff(req_cols, colnames(df_with_meta))) > 0) {
    stop(paste0(
      "there is a grouping condition column which is ",
      "not in the dataframe colnames"
    ))
  }

  req_grouping_conditions = c('CHR', 'POS', 'REF_Allele', 'ALT1')

  if (length(setdiff(req_grouping_conditions, grouping_conditions)) > 0){
    stop(paste0("This function currently requires at least the following ",
                "grouping conditions: ",
                paste(req_grouping_conditions, collapse = ","), ". ",
                "You may include additional, but these must be included."))
  }

  df_with_meta %>%
    # convert character vector to a vector of symbols,
    # then unquote that vector of symbols with !!!
    # https://stackoverflow.com/a/52572383
    dplyr::group_by(!!!dplyr::syms(grouping_conditions)) %>%
    dplyr::summarize(
      RealDepth = sum(RealDepth),
      Depth = sum(Depth),
      Reference = sum(Reference),
      Alternative1 = sum(Alternative1),
      .groups = "keep"
    ) %>%
    tidyr::unite("sample", !!!dplyr::syms(sample_columns), remove = FALSE) %>%
    dplyr::group_by(sample) %>%
    dplyr::arrange(CHR, POS, REF_Allele, ALT1) %>%
    dplyr::ungroup()
}
