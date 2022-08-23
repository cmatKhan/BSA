#TODO complete documentation

#' @title binner2 -- another binner
#' @description This is a function to make a new table with different
#'   bin sizes and the average of ∆% observed in the positions inside the bin.
#'   You can play with the bin size, and it uses the output of percenter
#'   as input.
#'
#' @note a bin is considered significant if half or more of the
#'   snps in the bin are significant
#'
#' @param input PARAM_DESCRIPTION
#' @param bin.size PARAM_DESCRIPTION, Default: 10000
#'
#' @return OUTPUT_DESCRIPTION
#'
#' @details DETAILS
#'
#' @export
binner2 <- function(input, bin.size = 10000) {

  # TODO replace with BSgenome
  KN99_gen <- data.frame(
    # note that the mito chrom is actually M in daniels data
    Chr = as.character(c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
                         "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                         "chr13", "chr14", "chrM")),
    size = c(2291500, 1621676, 1574972, 1084805, 1814975, 1422463, 1399209,
             1398693, 1186813, 1059962, 1562107, 774060, 756017, 942472, 24923)
  )
  input$CHROM <- as.character(input$CHROM)

  klaus <- data.frame(
    CHROM = as.character(), bin_floor = as.numeric(),
    bin_ceiling = as.numeric(), bin_middle = as.numeric(),
    deltaSNP = as.numeric(), triCubeDeltaSNP = as.numeric(),
    qvalue = as.numeric(), mut_perBin = as.numeric(), binSize = as.numeric(),
    Significance = as.logical()
  )
  for (i in 1:nrow(KN99_gen)) {
    bin.floor <- 0
    for (z in 1:ceiling(KN99_gen[i, 2] / bin.size)) {
      partial <- data.frame()
      bin.ceiling <- bin.floor + bin.size
      # check number of snps in bin
      partial <-subset(input, CHROM == KN99_gen$Chr[i] &
                         POS > bin.floor & POS <= bin.ceiling)
      # if 0, return NA
      if (nrow(partial) == 0) {
        deltaSNP <- NA
        triCubeDeltaSNP <- NA
      # if 1, return the sig of the snp
        } else if (nrow(partial) == 1) {
        deltaSNP <- partial$deltaSNP
        triCubeDeltaSNP <- partial$tricubeDeltaSNP
      # else return mean
        } else {
        deltaSNP <- mean(partial$deltaSNP, na.rm = T)
        triCubeDeltaSNP <- mean(partial$tricubeDeltaSNP, na.rm = T)
      }

      significance <- F
      # If half or more the snps in the window
      significance <-
        length(which(partial$qvalue < 0.1)) / length(partial$qvalue) >= 0.5
      if (is.na(significance)) {
        significance <- F
      }
      # significance=any(partial$qvalue<0.01, na.rm=T)
      # This is the case for window has one significant snp, all window is
      # significant.

      bin.middle <- mean(c(bin.floor, bin.ceiling))
      qvalue <- mean(partial$qvalue)
      klaus <- rbind(klaus, c(
        as.character(KN99_gen$Chr[i]),
        as.numeric(bin.floor),
        as.numeric(bin.ceiling),
        as.numeric(bin.middle),
        as.numeric(deltaSNP),
        as.numeric(triCubeDeltaSNP),
        as.numeric(qvalue),
        nrow(partial),
        bin.size,
        significance),
        stringsAsFactors = F)

      bin.floor <- bin.ceiling
    }
  }
  colnames(klaus) <- c("CHROM", "binFloor", "binCeiling", "binMiddle",
                       "deltaSNP", "smoothedDeltaSNP", "qvalue",
                       "mut_perBin", "Bin_size", "Significance")

  klaus[, 2:9] <- sapply(klaus[, 2:9], as.numeric)

  return(klaus)
}

