#' @title Another binner
#' @description see binner and binner2
#'
#' @param input PARAM_DESCRIPTION
#' @param bin.size PARAM_DESCRIPTION, Default: 10000
#'
#' @return OUTPUT_DESCRIPTION
#'
#' @details DETAILS
#'
#' @export
binner3 <- function(input, bin.size = 10000) {

  # TODO replace with BSgenome
  KN99_gen <- data.frame(
    Chr = c(
      "chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
      "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
      "chrM"
    ),
    size = c(
      2291500, 1621676, 1574972, 1084805, 1814975, 1422463, 1399209,
      1398693, 1186813, 1059962, 1562107, 774060, 756017, 942472, 24923
    )
  )

  smoothedDeltaSNP <- vector()
  for (i in 1:nrow(KN99_gen)) {
    bin.floor <- 0
    for (z in 1:ceiling(KN99_gen[i, 2] / bin.size)) {
      bin.ceiling <- bin.floor + bin.size
      partial <-
        subset(input,
               CHROM == KN99_gen$Chr[i] & POS > bin.floor & POS <= bin.ceiling)
      tempValue <- vector()
      if (nrow(partial) == 1) {
        tempValue <- partial$deltaSNP
      } else if (nrow(partial) > 1) {
        tempValue <- BSA::tricubeStat(POS = partial$POS,
                                 Stat = partial$deltaSNP,
                                 windowSize = bin.size)
      }
      if (length(tempValue) > 0) {
        smoothedDeltaSNP <- c(smoothedDeltaSNP, tempValue)
      }
      bin.floor <- bin.ceiling
    }
  }

  return(smoothedDeltaSNP)
}
