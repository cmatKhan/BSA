# TODO complete documentation

#' @title the tabler
#' @description FUNCTION_DESCRIPTION
#'
#' @param VCF PARAM_DESCRIPTION
#' @param i integer. Samples in a vcf file start in the 10th column.
#'   Tabler takes the information from the (9+i)th column from the vcf.
#'   Default: 1
#' @param minDepth Numeric. The minimal depth required at each
#'   position to call a different allele. Default: 5
#' @param cutoff Numeric. Minimum ratio of reads (as a number between 0 and 1)
#'   that supports either allele required to call an allele. Default: 0.75.
#' @param filtering Logical. If True, will fill up the column "filtGenotype".
#'   This is specially useful for single strains, not so much for BSA.
#'   The function is much faster if turned off. Default: F.
#' @param cutoff50 PARAM_DESCRIPTION, Default: 0.1
#'
#' @return OUTPUT_DESCRIPTION
#'
#' @details DETAILS
#'
#' @export
#'
#' @importFrom stringr str_split_fixed
#' @importFrom dplyr mutate
tabler <- function(VCF, i = 1, minDepth = 5, cutoff = 0.75,
                   filtering = F, cutoff50 = 0.1) {

  # This splits the FORMAT itens from the vcf into columns.
  c <- (stringr::str_split_fixed(VCF[, 9 + i], ":", n = 8))

  colnames(c) <- c("Genotype_GT", "Depth_DP", "AD", "Reference_RO", "QR", "Alternative_AO", "QA", "GL")

  # The two lines bellow convert any "." in the object
  # c(columns 3, 8 and the rest, respectively) to "0,0".
  c[, 3] <- sapply(c[, 3], USE.NAMES = F, FUN = function(x) {
    if (x == ".") {
      x <- "0,0"
    } else {
      x <- x
    }
  })
  c[, 8] <- sapply(c[, 8], USE.NAMES = F, FUN = function(x) {
    if (x == ".") {
      x <- "0,0"
    } else {
      x <- x
    }
  })
  c <- apply(c, c(1, 2), FUN = function(x) {
    if (x == ".") {
      x <- 0
    } else {
      x <- x
    }
  })

  # This gets the the number of Reference reads and ALternative
  # reads from the column AD from the object c.
  # This splits the two alternative allele counts into different columns.
  d <- stringr::str_split_fixed(c[, 3], ",", n = 3)

  # Convert c and d to numeric
  # Columns AD and GL will be converted to NA, generating errors.
  # But that's ok, we've got what we needed from them.
  c <- apply(c, c(1, 2), as.numeric)
  d <- apply(d, c(1, 2), as.numeric)

  # This gets the possible alternative alleles.
  # The split only happens when we have a second alternative allele.
  e <- stringr::str_split_fixed(VCF[, 5], ",", n = 3)

  data2 <- data.frame(
    CHR = VCF[, 1], POS = VCF[, 2], Ref_allele = VCF[, 4], ALT1 = e[, 1],
    ALT2 = e[, 2], ALT3 = e[, 3], Genotype = as.numeric(c[, 1]),
    Depth = as.numeric(c[, 2]), Reference = as.numeric(c[, 4]),
    Alternative1 = d[, 2], Alternative2 = d[, 3],
    filtGenotype = rep(".", nrow(e)), stringsAsFactors = FALSE
  )
  data2$Alternative2 <-
    sapply(data2$Alternative2, USE.NAMES = F, FUN = function(x) {
    if (is.na(x)) {
      x <- as.numeric(0)
    } else {
      x <- x
    }
  })

  data2$filtGenotype <- as.character(data2$filtGenotype)
  data2 <- dplyr::mutate(data2,
                         realDepth = Reference + Alternative1 + Alternative2,
                         Ref_percentage = Reference / realDepth,
                         Alt1_percentage = Alternative1 / realDepth,
                         Alt2_percentage = Alternative2 / realDepth
  )
  if (filtering == T) {
    for (j in 1:nrow(data2)) {
      # Loop for defining the filtered Genotype. Filter by minimum depth and
      # Genotype cutoff. The cutoff is the minimum percentage of reads
      # necessary to define a genotype. Default = 0.75. the cutoff50 is the
      # distance to 0.5 to a position so it can be called as
      # "Mixed", as in mixed strains.
      data2$filtGenotype
      if (is.na(data2$Reference[j])) {
        x <- "noReads"
      } else if (data2$realDepth[j] < (minDepth + 1)) {
        data2[j, 12] <- "lowDepth"
      } else if (data2$Reference[j] > cutoff * data2$realDepth[j]) {
        data2[j, 12] <- "Reference"
      } else if (data2$Alternative1[j] > cutoff * data2$realDepth[j]) {
        data2[j, 12] <- "Alternative1"
      } else if (data2$Alternative2[j] > cutoff * data2$realDepth[j]) {
        data2[j, 12] <- "Alternative2"
      } else if (
        (data2$Ref_percentage[j] > (0.5 - cutoff50) &
         data2$Ref_percentage[j] < (0.5 + cutoff50)) &
        (data2$Alt1_percentage[j] > (0.5 - cutoff50) &
         data2$Alt1_percentage[j] < (0.5 + cutoff50))) {
        data2[j, 12] <- "Mixed"
      } else if ((data2$Alternative2[j] <= cutoff * data2$realDepth[j] &
                  data2$Alternative1[j] <= cutoff * data2$realDepth[j] &
                  data2$Reference[j] <= cutoff * data2$realDepth[j])) {
        data2[j, 12] <- "undefinedGenotype"
      }
    }
  }

  rm(c)
  rm(d)

  return(data2)
}


