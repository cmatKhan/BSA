#' @param ... Other arguments passed to \code{\link[locfit]{locfit}} and
#' @param altFreq1 numeric. The alternate allele frequency for bulk A.
#' @param altFreq2 numeric. The alternate allele frequency for bulk B.
#' @param bin.size The size of the window to be used in the smoothing.
#' @param bulkSize
#' @param cutoff Numeric. Minimum ratio of reads (as a number between 0 and 1) that supports either allele required to call an allele. Default=0.75.
#' @param deltaSNP a vector of delta SNP values for use for QTL region filtering
#' @param depth integer. A read depth for which to replicate SNP-index calls.
#' @param filter numeric. an optional minimum SNP-index filter
#' @param filtering Logical. If True, will fill up the column "filtGenotype". This is specially useful for single strains, not so much for BSA. The function is much faster if turned off. Default=F.
#' @param filterThreshold The absolute delta SNP index to use to filter out
#' @param filterThreshold The absolute delta SNP index to use to filter out putative QTL
#' @param Gprime a vector of G prime values (tricube weighted G statistics)
#' @param HighAlt A vector of the alternate allele depth in the high bulk
#' @param highBulk Character. The highbulk of the comparison made.
#' @param HighRef A vector of the reference allele depth in the high bulk
#' @param i integer. Samples in a vcf file start in the 10th column. Tabler takes the information from the (9+i)th column from the vcf.
#' @param input
#' @param LowAlt A vector of the alternate allele depth in the low bulk
#' @param lowBulk Character. The lowbulk of the comparison made.
#' @param LowRef A vector of the reference allele depth in the low bulk
#' @param Max depth: highr fith percentile > quantile(sum of depth in both bulks, na.rm = T, probs = 0.95))
#' @param maxDepthPercentile depth: higher fifth percentile > quantile(sum of depth in both bulks, na.rm = T, probs = 0.95))
#' @param Min sample depth - Half the min depth
#' @param minDepth Numeric. The minimal depth required at each position to call a different allele. Default=5.
#' @param minDepthPercentile sample depth - Half the min depth
#' @param MinDepthPercentile: lower tenth percentile > quantile(sum of depth in both bulks, na.rm = T, probs = 0.1))
#' @param minimumDepthPerAllele
#' @param outlierFilter one of either "deltaSNP" or "Hampel". Method for
#' @param POS A vector of genomic positions for each SNP
#' @param replicates integer. The number of bootstrap replications.
#' @param SNPcomparison: lower tenth percentile > quantile(sum of depth in both bulks, na.rm = T, probs = 0.1))
#' @param SNPset  List. The main input for this function is a list in which all
#' @param SNPset Data frame SNP set containing previously filtered SNPs
#' @param Stat A vector of values for a given statistic for each SNP
#' @param VCF
#' @param windowSize
#' @param windowSize the window size (in base pairs) bracketing each SNP for which
#' @param WinSize the window size (in base pairs) bracketing each SNP for which
