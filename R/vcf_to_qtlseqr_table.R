# TODO INCLUDE PARENT FILTER

#' @title Set Filt_Genotype column
#' @description internal helper function to create Filt_Genotype column. note
#'   this column name is from Daniel's original BSA analysis. It just means
#'   a column which calls genotype based on frequency of allele with RealDepth,
#'   which are only those bases which pass the variant caller thresholds and
#'   exist (not missing in the CIGAR), in the denominator
#'
#' @param depth depth of a given variant, at a given position
#' @param ref_freq frequency of a given variant, at a given position
#' @param depth_thres if less than depth threshold, or is.na(depth), return
#'   value is lowDepth
#' @param ref_freq_thres if ref is greater than or equal to this number, and
#'   depth passes threshold, call genotype Reference. note that genotype is called
#'   Alternative if passes depth threshold and ref_freq is less than or equal to
#'   1 - ref_freq_thres
#'
#' @return one of undefinedGenotype, lowDepth, Reference, Alternative
set_genotype = function(depth, ref_freq, depth_thres, ref_freq_thres){


  genotype = 'undefinedGenotype'

  if(is.na(depth) | depth < depth_thres){
    genotype = 'lowDepth'
  } else if(ref_freq >= ref_freq_thres){
    genotype = 'Reference'
  } else if(ref_freq <= 1-ref_freq_thres){
    genotype = 'Alternative'
  }
  genotype
}

#' @title Tidy data from gds file
#' @description Tidy data to long format with column sample, variant, and the
#'   metric name. internal use only
#'
#' @param data_mat a data matrix in sample x variant format
#' @param values_cname name of the metric -- eg RealDepth
#' @param rnames vector to name the columns, eg sample.ids
#' @param cnames vector to name the rows, eg variant_ids
#'
#' @return a tidy dataframe in long format with columns sample, variant, metric
#'
#' @importFrom dplyr as_tibble
#' @importFrom tidyr pivot_longer
tidy_metrics = function(data_mat, values_cname, rnames, cnames){

  errors = list(
    dimension_error = paste0("vcf_to_qtlseqr_table error: rownames do not ",
                             "match data matrix. see tidy_metrics() function")
  )

  if(length(rnames) != dim(data_mat)[1] | length(cnames) != dim(data_mat)[2]){
    stop(errors$dimension_error)
  }

  data_mat %>%
    `rownames<-`(rnames) %>%
    `colnames<-`(cnames) %>%
    dplyr::as_tibble(rownames = "sample") %>%
    tidyr::pivot_longer(-sample, names_to = "variant", values_to = values_cname)
}


#' @title VCF to QTLseqR table
#' @description Process the VCF, create a gds file, and return a table
#' which is in long format. See return for column list
#'
#' @param vcf_path Path to a VCF file, presumably one with a number of samples
#' @param gds_outdir directory to which to write the gds file
#' @param depth_thres minimum required (filtered) depth to consider calling
#'   a genotype. Less than this and the genotype is labelled lowDepth, Default: 5
#' @param ref_freq_thres minimum (greater than or equal to) required percentage
#'   for a genotype to be called Reference or Alternative, Default: 0.9
#' @param parent_filter Boolean. Set to TRUE to retain only those loci which
#'   overlap parent strains, KN99alpha and TDY1993
#' @param single_allele_loci_only Boolean, set to TRUE to exclude all
#'   multi allelic loci. Set to FALSE to keep all variants, Default: TRUE
#' @param overwrite_gds if the gds already exists, skip creating and just open.
#'   Default is FALSE
#' @param verbose Boolean. Set to true to set the SeqArray functions to verbose.
#'   Default is FALSE.
#'
#' @return A dataframe with the following columns:  CHR POS REF_Allele ALT1 QUAL
#'   Depth variant sample RealDepth Reference Alternative1 genotype
#'   Ref_percentage Alt1_percentage Filt_Genotype
#'
#' @details This function is a re-interpretation of Daniel's script
#'   /scratch/mblab/daniel.agustinho/tools/VCF_tabler.sh.
#'   This script parsed the VCF files through a series of awk commands.
#'   Critically, the info line for the VCF is expected to be in the following
#'   format: GT:DP:AD:RO:QR:AO:QA:GL. If you read that script, you'll see
#'   that the awk command extracted columns 1-9, which are
#'   CHROM  POS ID  REF ALT QUAL  FILTER  INFO  FORMAT
#'   it a loop, it also extracts a given sample column. That column format is
#'   given by the INFO column, and is GT:DP:AD:RO:QR:AO:QA:GL, as stated above.
#'   Real depth was calculated by replacing the colons with field separators
#'   and extracting columns 13 and 15, which correspond to RO and AO respectively.
#'   This is mimicked in this re-interpreted function.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[SeqArray]{seqVCF2GDS}},
#'  \code{\link[SeqArray]{seqGetData}}
#'
#' @rdname vcf_to_qtlseqr_table
#'
#' @export
#' @importFrom stringr str_remove_all
#' @importFrom SeqArray seqVCF2GDS seqOpen seqGetData seqSetFilter seqResetFilter seqClose
#' @importFrom SeqVarTools nAlleles
#' @importFrom purrr map2 reduce
#' @importFrom dplyr left_join tibble mutate arrange
vcf_to_qtlseqr_table = function(vcf_path, gds_outdir,
                                depth_thres=5, ref_freq_thres=.9,
                                single_allele_loci_only = TRUE,
                                parent_filter = TRUE,
                                overwrite = FALSE,
                                verbose = FALSE){

  # construct output path for gds file
  vcf_name = stringr::str_remove_all(basename(vcf_path), ".vcf$|.vcf.gz")
  gds_path = file.path(gds_outdir, paste0(vcf_name, ".gds"))

  # if overwrite is false and the gds already exists, skip creating the gds
  if(!file.exists(gds_path) | overwrite == TRUE){
    message("Creating GDS from VCF...")
    # create GDS file from VCF
    SeqArray::seqVCF2GDS(vcf_path,out.fn = gds_path, verbose = verbose)

  }

  # open the gds file
  message("opening gds...")
  gds = SeqArray::seqOpen(gds_path)

  # this is wrapped in a tryCatch block in order to ensure that the gds object
  # is closed, even if there is an error
  tryCatch({
    # if single_allele_loci_only, set a filter on the gds object such that only
    # single allele loci are extracted in following steps
    if(single_allele_loci_only){
      message("Setting single allele loci filter...")
      single_allele_fltr = SeqArray::seqGetData(gds, 'variant.id')[SeqVarTools::nAlleles(gds) == 2]
      message(sprintf("...total variants: %s", length(SeqVarTools::nAlleles(gds))))
      message("...FALSE refers to the number of multi allelic loci")
      print(table(SeqVarTools::nAlleles(gds) == 2))
      SeqArray::seqSetFilter(gds, variant.id=single_allele_fltr)
    }

    # slots from which to extract data from the gds object. See seqSummary(gds)
    # for details on what columns we may get. Also see
    # https://gatk.broadinstitute.org/hc/en-us/articles/360036711711-DepthPerAlleleBySample
    # and
    # https://gatk.broadinstitute.org/hc/en-us/articles/360036347832-Coverage
    message("Extracting data from gds...")
    data_to_extract = c("chromosome", "position", "$ref", "$alt",
                        "annotation/qual", "genotype", "annotation/info/DP",
                        "annotation/format/DP", "annotation/format/RO",
                        "annotation/format/AO", "sample.id")
    # extract the data
    gds_extracted_data = SeqArray::seqGetData(gds, data_to_extract)

    # note: this must extract data over the 'unfiltered' depth -- values are
    # very, very different from the "realdepth" calculation below
    # ref_allele_freq = seqAlleleFreq(gds)
  }, error = function(e){
    stop("Error extracting data from gds object: ", e)
  # whether or not there is an error, close the gds object
  }, finally = {
    SeqArray::seqClose(gds)
  })

  # create unique IDs for the variants
  variant_ids = paste0("var_", seq(1,dim(gds_extracted_data$`annotation/format/DP`)[2]))

  # these metrics are in sample x variant matricies. The following two steps
  # transform them into tidy, long format dataframes
  message("Transforming depth data...")
  data_to_transform = list(
    genotype     = gds_extracted_data$genotype[1,,],
    RealDepth    = gds_extracted_data$`annotation/format/DP`,
    Reference    = gds_extracted_data$`annotation/format/RO`,
    Alternative1 = gds_extracted_data$`annotation/format/AO`$data)
  # result is a long data frame with columns sample, variant, RealDepth,
  # Reference (which is the reference allele depth) and Alternative1, which is
  # the Alternative depth. Note that if single_allele_loci_only is not set,
  # then this is not actually Alternative1, but just Alternative
  depth_metrics_df = suppressMessages(purrr::map2(data_to_transform,
                                 names(data_to_transform),
                                 tidy_metrics,
                                 gds_extracted_data$sample.id, variant_ids) %>%
    purrr::reduce(dplyr::left_join))

  # create a vector which represents chr_pos of variants in
  # the parent strains
  if(parent_filter){
    # CREATE PARENT FILTER
    # ```{r, eval = FALSE}
    #
    # kn99a_tdy1993_variant_filter = parents$KN99a$Filt_Genotype == "Reference" &
    #   parents$TDY1993$Filt_Genotype == "Alternative"
    # parent_variants = parents$KN99a[kn99a_tdy1993_variant_filter,]
    # ```
    #
    # # Gather variant strain (pool) data
    #
    # `mySamples` will be a list of dataframes. Each frame is (row) filtered for
    # only those positions contained in `parent_variants`
    #
    # ```{r, eval = FALSE}
    # poolsBSA2 <- list.files(
    #   path = "/mnt/lts/personal/chasem/DNA/EXP#031/VCF_tables") %>%
    #   sub(pattern = ".txt", replacement = "")
    #
    #
    # samples_df_bsa2 = foreach(
    #   i = poolsBSA2,
    #   .inorder = TRUE,
    #   .combine = 'rbind'
    # ) %do% {
    #   df = read.table(file = paste0(
    #     "/mnt/lts/personal/chasem/DNA/EXP#031/VCF_tables/", i, ".txt"),
    #     header = T,
    #     stringsAsFactors = F,
    #     numerals = "allow.loss")
    #
    #   # only retain those variants which overlap kn99a vs tdy1993 variants
    #   parent_filter = paste(df$CHR, df$POS, sep = "_") %in%
    #     paste(parent_variants$CHR, parent_variants$POS, sep = "_")
    #
    #   df = df[parent_filter, ]
    #
    #   df$Alternative1 <- as.integer(df$Alternative1)
    #
    #   df$sample = i
    #
    #   df
    # }
    #
    warning(paste0("PARENT FILTER NOT IMPLEMENTED -- ",
                   "SETTING TO FALSE. NO POSITIONS ",
                   "WILL BE FILTERED.\n"))
    parent_filter = FALSE
    parent_filter_vector = NULL
  }

  # reshape all of the extracted data into a tidy, long data frame.
  # note that columns are named and formatted according to Daniel's
  # original table specification
  message("Creating QTLseqR table...")
  suppressMessages(dplyr::tibble(CHR = gds_extracted_data$chromosome,
                POS = gds_extracted_data$position,
                REF_Allele = gds_extracted_data$`$ref`,
                ALT1 = gds_extracted_data$`$alt`,
                QUAL = gds_extracted_data$`annotation/qual`,
                Depth = gds_extracted_data$`annotation/info/DP`,
              # see the note above on ref_allele_freq -- this and the RealDepth
              # do not agree
              # Ref_percentage = ref_allele_freq,
              # Alt1_percentage = 1-ref_allele_freq
  ) %>%
    dplyr::mutate(variant = variant_ids) %>%
    dplyr::left_join(depth_metrics_df) %>%
    dplyr::mutate(Ref_percentage = Reference/RealDepth) %>%
    dplyr::mutate(Alt1_percentage = 1-Ref_percentage) %>%
    dplyr::mutate(Filt_Genotype = unlist(purrr::map2(RealDepth, Ref_percentage,
                                       set_genotype,
                                       depth_thres, ref_freq_thres ))) %>%
    dplyr::arrange(sample)) %>%
    # the following two lines apply the parent filter if it is set to TRUE
    # otherwise, the dataframe is unfiltered at this point
    unite("cp", CHR, POS, sep ="_", remove=FALSE) %>%
    filter(if(parent_filter) cp %in% parent_filter_vector else TRUE) %>%
    dplyr::select(-cp)
}
