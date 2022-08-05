# TODO complete documentation

#' @title a percenter function
#' @description this function will make a new table with the differences
#'   in percentages (between YPD - inoculum, and Lung - inoculum)
#'   in a sample using the tabler output list as input. The sampleColumn
#'   should be the column number that shows the sample's inoculum.

#'
#' @param input PARAM_DESCRIPTION, Default: allSamples
#' @param parameter PARAM_DESCRIPTION, Default: 'Reference'
#' @param sampleColumn PARAM_DESCRIPTION, Default: 1
#'
#' @return OUTPUT_DESCRIPTION
#'
#' @details DETAILS
#'
#' @export
#'
#' @importFrom generics as.factor
percenter <- function(input = allSamples, parameter = "Reference", sampleColumn = 1) {
  # this function makes a list with the differences in percentage between the Lung/YPD samples and the inoculum. The parameter for testing is
  # either the reference or the alternative alleles.
  v <- sampleColumn
  reg <- input[[v]][, c(1:6)]
  if (parameter == "Reference") {
    reg$Lung <- input[[v + 1]]$Ref_percentage - input[[v]]$Ref_percentage
    reg$YPD <- input[[v + 2]]$Ref_percentage - input[[v]]$Ref_percentage
  } else
  if (parameter == "Alternative1") {
    reg$Lung <- input[[v + 1]]$Alt1_percentage - input[[v]]$Alt1_percentage
    reg$YPD <- input[[v + 2]]$Alt1_percentage - input[[v]]$Alt1_percentage
  } else
  if (parameter == "Alternative2") {
    reg$Lung <- input[[v + 1]]$Alt2_percentage - input[[v]]$Alt2_percentage
    reg$YPD <- input[[v + 2]]$Alt1_percentage - input[[v]]$Alt1_percentage
  }
  reg$CHR <- generics::as.factor(reg$CHR)
  return(reg)
}
