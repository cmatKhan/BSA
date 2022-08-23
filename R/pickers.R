# TODO create a class which instantiates with the sample_metadata and vcf_df,
#   possesses a method to extract groups

#' @title Transform VCF Table for QTLseqR
#' @description Function to set comparisons between two samples and output
#'   the table in a way that the QTLseq package can use. Uses two element of
#'   a list to create a table that can be used by the package QTLseq.
#'
#' @param samples_with_meta_df A dataframe representing samples
#'   \code{\link[BSA]{vcf_to_qtlseqr}} joined with the sample metadata
#' @param low_bulk_sample String. The name of the sample (and of the element in
#'   the vcf_df list) corresponding to the Low Bulk.
#' @param high_bulk_sample The name of the sample (and of the element in
#'   the vcf_df list) corresponding to the High Bulk.
#' @param sample_col name of the column with the sample names. Must contain
#'   the low_bulk_sample and high_bulk_sample. Default 'sample'.
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
picker2 = function(samples_with_meta_df,
                    low_bulk_sample, high_bulk_sample,
                    sample_col = "sample",
                    bulk_col = "bulk") {

  if(length(setdiff(c(sample_col, bulk_col),
                    colnames(samples_with_meta_df))) > 0){
    stop(paste0(sprintf("sample_col: %s and bulk_col: %s ", sample_col, bulk_col ),
                 "must be in the colnames of the input dataframe."))
  }

  if (!low_bulk_sample %in% samples_with_meta_df[[sample_col]] |
      !high_bulk_sample %in% samples_with_meta_df[[sample_col]]) {
    stop(paste0("Confirm if the bulks' names are present",
         "in the SNP set (vcf file) provided."))
  }

  # after joining, reduce the joined dataframe down to just these columns
  select_cols = c('CHR', 'POS', 'REF_Allele', 'ALT1',
                  'Depth', 'RealDepth', 'Reference', 'Alternative1',
                  bulk_col)

  if(length(setdiff(select_cols,
                    colnames(samples_with_meta_df))) > 0){
    stop(paste0("the following columns must be in the input dataframe: ",
                paste(select_cols, collapse=",")))
  }

  # the column on which to pivot
  names_from_col = "bulk"

  # column on which to filter the low_bulk and high_bulk samples
  sample_name_col = "sample"

  # it is from these columns which the values -- entries -- of the new "wider"
  # columns are created. If the levels of the `names_from_col` are {low, high},
  # then the new columns will be Genotype_low, Depth_low, ..., Genotype_high,
  # Depth_high, ...
  # Note that another way of doing this might be to use the - operator to
  # exclude columns from the pivot -- haven't tried this.
  value_cols = c('Depth', 'RealDepth', 'Reference', 'Alternative1')
  if(length(setdiff(value_cols, select_cols)) > 0){
    stop("See picker2 code: value_cols must be a subset of select_cols")
  }

  # rename columns from old to new. Note that the index of the column in old
  # must correspond to the index of the column in new
  # NOTE that the `mutate` function below depends on these names -- if this
  # is changed, then the corresponding field in `mutate` needs to be addressed
  rename_cols = list(

    old = c("CHR","REF_Allele","ALT1",
            "Reference_low","Alternative1_low","RealDepth_low",
            "Reference_high","Alternative1_high","RealDepth_high"),

    new = c("CHROM","REF","ALT",
            "AD_REF.LOW","AD_ALT.LOW", "DP.LOW",
            "AD_REF.HIGH","AD_ALT.HIGH","DP.HIGH")
  )

  # transform the dataframe
  samples_with_meta_df %>%
    filter(!!rlang::sym(sample_name_col) %in%
             c(low_bulk_sample, high_bulk_sample)) %>%
    select(all_of(select_cols)) %>%
    pivot_wider(names_from = !!rlang::sym(names_from_col),
                values_from = all_of(value_cols)) %>%
    rename_with(~rename_cols$new, all_of(rename_cols$old)) %>%
    # add columns -- NOTE that these are hard coded. If rename_cols is changed,
    # this must be changed
    mutate(GQ.LOW             = NA,
           PL.LOW             = NA,
           SNPindex.LOW       = AD_ALT.LOW / DP.LOW,
           GQ.HIGH             = NA,
           PL.HIGH             = NA,
           SNPindex.HIGH       = AD_ALT.HIGH / DP.HIGH,
           REF_FRQ            = (AD_ALT.LOW + AD_REF.HIGH) / (DP.LOW + DP.HIGH),
           deltaSNP           = SNPindex.HIGH - SNPindex.LOW
    )

}
