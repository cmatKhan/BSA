#TODO better title. Review docs and improve

#' @title pick something
#' @description Function to set comparisons between two samples and output
#'   the table in a way that the QTLseq package can use. Uses two element of
#'   a list to create a table that can be used by the package QTLseq.
#'
#' @param SNPset  List. Each eelement of the list is a table showing data
#'   from one sample, the table is set as the output of the tabler function.
#' @param lowBulk String. The name of the sample (and of the element in
#'   the SNPset list) corresponding to the Low Bulk.
#' @param highBulk The name of the sample (and of the element in
#'   the SNPset list) corresponding to the High Bulk.
#'
#' @return a data.frame with the format that can be understood by QTLseq.
#'
#' @details This will be in the format used
#'   by QTLseq package. Picker2 has some adaptations from picker 1 for
#'   the tabler function in awk. RealDepth instead of realDepth,
#'   and no Alt2 allele.
#'
#' @export
#'
#' @importFrom dplyr mutate
picker2 <- function(SNPset, lowBulk, highBulk) {

  if (!lowBulk %in% names(SNPset) | !highBulk %in% names(SNPset)) {
    stop(paste0("Confirm if the bulks' names are present",
         "in the SNP set (vcf file) provided."))
  }

  metaTable <- data.frame(
    CHROM = SNPset[[lowBulk]]$CHR,
    POS = SNPset[[lowBulk]]$POS,
    REF = SNPset[[lowBulk]]$Ref_Allele,
    ALT = SNPset[[lowBulk]]$ALT1,
    AD_REF.LOW = SNPset[[lowBulk]]$Reference,
    AD_ALT.LOW = (SNPset[[lowBulk]]$Alternative1),
    DP.LOW = SNPset[[lowBulk]]$RealDepth,
    GQ.LOW = rep(NA, nrow(SNPset[[lowBulk]])),
    PL.LOW = rep(NA, nrow(SNPset[[lowBulk]])),
    SNPindex.LOW = (as.numeric(SNPset[[lowBulk]]$Alternative1)) /
      as.numeric(SNPset[[lowBulk]]$RealDepth),
    AD_REF.HIGH = SNPset[[highBulk]]$Reference,
    AD_ALT.HIGH = (SNPset[[highBulk]]$Alternative1),
    DP.HIGH = SNPset[[highBulk]]$RealDepth,
    GQ.HIGH = rep(NA, nrow(SNPset[[highBulk]])),
    PL.HIGH = rep(NA, nrow(SNPset[[highBulk]])),
    SNPindex.HIGH = (as.numeric(SNPset[[highBulk]]$Alternative1)) /
      as.numeric(SNPset[[highBulk]]$RealDepth), # AD_ALT.HIGH/DP.HIGH,
    REF_FRQ = ((SNPset[[highBulk]]$Reference + SNPset[[lowBulk]]$Reference) /
                 (SNPset[[highBulk]]$RealDepth + SNPset[[lowBulk]]$RealDepth))
  )
  metaTable <- dplyr::mutate(metaTable, deltaSNP = SNPindex.HIGH - SNPindex.LOW)
  return(metaTable)
}
