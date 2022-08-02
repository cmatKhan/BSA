#' @title pick something
#' @description FUNCTION_DESCRIPTION
#'
#' @param SNPset PARAM_DESCRIPTION
#' @param lowBulk PARAM_DESCRIPTION
#' @param highBulk PARAM_DESCRIPTION
#'
#' @return OUTPUT_DESCRIPTION
#'
#' @details DETAILS
#'
#' @export
#'
#' @importFrom dplyr mutate
picker2 <- function(SNPset, lowBulk, highBulk) {

  if (!lowBulk %in% names(SNPset) | !highBulk %in% names(SNPset)) {
    stop("Confirm if the bulks' names are present in the SNP set (vcf file) provided.")
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
    SNPindex.LOW = (as.numeric(SNPset[[lowBulk]]$Alternative1)) / as.numeric(SNPset[[lowBulk]]$RealDepth),
    AD_REF.HIGH = SNPset[[highBulk]]$Reference,
    AD_ALT.HIGH = (SNPset[[highBulk]]$Alternative1),
    DP.HIGH = SNPset[[highBulk]]$RealDepth,
    GQ.HIGH = rep(NA, nrow(SNPset[[highBulk]])),
    PL.HIGH = rep(NA, nrow(SNPset[[highBulk]])),
    SNPindex.HIGH = (as.numeric(SNPset[[highBulk]]$Alternative1)) / as.numeric(SNPset[[highBulk]]$RealDepth), # AD_ALT.HIGH/DP.HIGH,
    REF_FRQ = ((SNPset[[highBulk]]$Reference + SNPset[[lowBulk]]$Reference) / (SNPset[[highBulk]]$RealDepth + SNPset[[lowBulk]]$RealDepth))
  )
  metaTable <- dplyr::mutate(metaTable, deltaSNP = SNPindex.HIGH - SNPindex.LOW)
  return(metaTable)
}
