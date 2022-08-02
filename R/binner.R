#
#' @title Binner -- the first
#' @description This is a function to make a new table with different bin sizes
#'   and the average of ∆% observed in the positions inside the bin.
#'   You can play with the bin size, and it uses the output of percenter as input.
#' @param input PARAM_DESCRIPTION
#' @param bin.size PARAM_DESCRIPTION, Default: 10000
#'
#' @return OUTPUT_DESCRIPTION
#'
#' @details DETAILS
#'
#' @export
#'
#' @importFrom timeSeries colnames
binner <- function(input, bin.size = 10000) {

  # TODO replace with BSgenome
  KN99_gen <- data.frame(
    Chr = c(
      "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
      "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chrM"
    ),
    size = c(
      2291500, 1621676, 1574972, 1084805, 1814975, 1422463, 1399209,
      1398693, 1186813, 1059962, 1562107, 774060, 756017, 942472, 24923
    )
  )

  klaus <- data.frame(
    CHR = as.character(), bin_floor = as.numeric(),
    bin_ceiling = as.numeric(), bin_middle = as.numeric(),
    lungAverage = as.numeric(), YPDaverage = as.numeric(),
    mut_perBin = as.numeric(), binSize = as.numeric()
  )
  for (i in 1:nrow(KN99_gen)) {
    bin.floor <- 0
    for (z in 1:ceiling(KN99_gen[i, 2] / bin.size)) {
      bin.ceiling <- bin.floor + bin.size
      partial <-
        subset(input,
               CHR == KN99_gen$Chr[i] & POS > bin.floor & POS <= bin.ceiling)
      lungAve <- mean(partial$Lung, na.rm = T)
      YPDAve <- mean(partial$YPD, na.rm = T)
      bin.middle <- mean(c(bin.floor, bin.ceiling))
      klaus <- rbind(klaus, c(
        as.character(KN99_gen$Chr[i]),
        as.numeric(bin.floor),
        as.numeric(bin.ceiling),
        as.numeric(bin.middle),
        as.numeric(lungAve),
        as.numeric(YPDAve),
        nrow(partial),
        bin.size),
        stringsAsFactors = F)

      bin.floor <- bin.ceiling
    }
  }

  timeSeries::colnames(klaus) <- c("Chr", "binFloor", "binCeiling", "binMiddle",
                       "lungMean", "YPDmean", "mut_perBin", "Bin_size")
  klaus[, 2:ncol(klaus)] <- sapply(klaus[, 2:ncol(klaus)], as.numeric)

  return(klaus)
}
