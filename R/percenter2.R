#' @title Another percenter function
#' @description FUNCTION_DESCRIPTION
#'
#' @param input PARAM_DESCRIPTION
#' @param parameter PARAM_DESCRIPTION, Default: 'Reference'
#' @param sampleColumn PARAM_DESCRIPTION, Default: 1
#'
#' @return OUTPUT_DESCRIPTION
#'
#' @details DETAILS
#'
#' @export
#'
#' @importFrom stringr str_split_fixed
#' @importFrom generics as.factor
percenter2 <- function(input, parameter = "Reference", sampleColumn = 1) {
  # this function makes a list with the differences in percentage between the Lung/YPD samples and the inoculum. The parameter for testing is
  # either the reference or the alternative alleles.
  v <- sampleColumn
  tester1 <- BSA::tabler(i = v)
  tester2 <- BSA::tabler(i = (v + 1))
  tester3 <- BSA::tabler(i = (v + 2))
  reg <- tester1[, c(1:6)]
  if (parameter == "Reference") {
    reg$Lung <- tester2$Ref_percentage - tester1$Ref_percentage
    reg$YPD <- tester2$Ref_percentage - tester1$Ref_percentage
  } else
  if (parameter == "Alternative1") {
    reg$Lung <- tester2$Alt1_percentage - tester1$Alt1_percentage
    reg$YPD <- tester2$Alt1_percentage - tester1$Alt1_percentage
  } else
  if (parameter == "Alternative2") {
    reg$Lung <- tester2$Alt2_percentage - tester1$Alt2_percentage
    reg$YPD <- tester2$Alt1_percentage - tester1$Alt1_percentage
  }

  sampleName <- colnames(tab)
  amp <- stringr::str_split_fixed(sampleName, ".I", n = 2)[9 + v]

  reg$replicate <- rep(amp, nrow(reg))
  rm(tester1)
  rm(tester2)
  rm(tester3)
  reg$CHR <- generics::as.factor(reg$CHR)
  reg$replicate <- generics::as.factor(reg$replicate)
  return(reg)
}
