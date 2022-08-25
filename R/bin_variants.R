
#' @title Bin Variants and calculate stats over bins
#' @description FUNCTION_DESCRIPTION
#'
#' @param variant_metrics_df PARAM_DESCRIPTION
#' @param tiled_genome_df PARAM_DESCRIPTION
#' @param chr_seqlengths PARAM_DESCRIPTION
#' @param chr_colname PARAM_DESCRIPTION
#' @param qvalue_lower_thres PARAM_DESCRIPTION, Default: 0.1
#' @param qvalue_upper_thres PARAM_DESCRIPTION, Default: 0.5
#'
#' @return OUTPUT_DESCRIPTION
#'
#' @details DETAILS
#'
#' @examples
#'
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{group_by}},
#'  \code{\link[dplyr]{group_split}},
#'  \code{\link[dplyr]{pull}},
#'  \code{\link[dplyr]{summarise}},
#'  \code{\link[dplyr]{context}},
#'  \code{\link[dplyr]{mutate}},
#'  \code{\link[dplyr]{select}},
#'  \code{\link[dplyr]{reexports}}
#'  \code{\link[rlang]{sym}}
#'  \code{\link[purrr]{map}}
#'  \code{\link[stringr]{str_remove}}
#'  \code{\link[tidyr]{separate}}
#'
#' @rdname bin_variants
#'
#' @export
#'
#' @importFrom dplyr group_by group_split pull summarize n mutate select all_of
#' @importFrom rlang sym
#' @importFrom purrr map
#' @importFrom stringr str_remove_all
#' @importFrom tidyr separate
bin_variants = function(variant_metrics_df,
                        tiled_genome_df,
                        chr_seqlengths,
                        chr_colname,
                        qvalue_lower_thres = 0.1,
                        qvalue_upper_thres = 0.5){

  if(!chr_colname %in% colnames(variant_metrics_df)){
    stop(sprintf("%s not in variant_metrics_df colnames", chr_colname))
  } else if(!chr_colname %in% colnames(tiled_genome_df)){
    stop(sprintf("%s not in tiled_genome_df colnames", chr_colname))
  }

  final_col_order = c(chr_colname, 'tile', 'binFloor', 'binCeiling',
                      'binMiddle', 'mean_deltaSNP', 'mean_smoothed_delta_snp',
                      'mean_qvalue', 'mut_perBin', 'bin_size', 'significance')

  variant_metrics_chr_split = variant_metrics_df %>%
    group_by(!!rlang::sym(chr_colname)) %>%
    group_split()

  names(variant_metrics_chr_split) = unlist(
    map(variant_metrics_chr_split, ~unique(pull(.,chr_colname)))
  )

  tiled_genome_chr_split = tiled_genome_df %>%
    group_by(!!rlang::sym(chr_colname)) %>%
    group_split()

  names(tiled_genome_chr_split) = unlist(
    map(tiled_genome_chr_split, ~unique(pull(.,chr_colname)))
  )

  out = map(names(variant_metrics_chr_split),
            ~tile_metrics(variant_metrics_chr_split[[.]],
                          tiled_genome_chr_split[[.]],
                          chr_seqlengths[[.]])) %>%
    do.call('rbind', .) %>%
    group_by(!!rlang::sym(chr_colname), tile) %>%
    summarize(mean_deltaSNP             = mean(deltaSNP, na.rm=TRUE),
              mean_smoothed_delta_snp   = mean(tricubeDeltaSNP, na.rm = TRUE),
              mean_qvalue               = mean(qvalue, na.rm = TRUE),
              qval_below_thres_n        = sum(qvalue < qvalue_lower_thres),
              qval_above_thres_n        = sum(qvalue >= qvalue_upper_thres),
              mut_perBin                = n(),
              .groups = 'keep') %>%
    mutate(significance = qval_below_thres_n / qval_above_thres_n,
           tile = str_remove_all(tile, '\\[|\\]|\\(|\\)')) %>%
    separate(tile, c('binFloor', 'binCeiling'), sep=",", remove = FALSE) %>%
    mutate(binFloor = as.integer(binFloor),
           binCeiling = as.integer(binCeiling)) %>%
    mutate(bin_size = binCeiling - binFloor,
           binMiddle = (binFloor + binCeiling) / 2)

  if(length(setdiff(final_col_order, colnames(out))) > 0){
    stop(paste0(
      sprintf("The following columns are expected in the output: %s.\n",
              paste(final_col_order, collapse = ",")),
      sprintf("But the output columns do not possess one or more: %s",
              paste(colnames(out), collapse=","))) )
  } else{
    out %>%
      dplyr:: select(dplyr::all_of(final_col_order))
  }

}

#' @title Tile Metrics Dataframe
#' @description FUNCTION_DESCRIPTION
#'
#' @param metrics_df PARAM_DESCRIPTION
#' @param genome_tile_df PARAM_DESCRIPTION
#' @param chr_seqlength PARAM_DESCRIPTION
#' @param cut_col PARAM_DESCRIPTION, Default: 'POS'
#' @param left_inclusive PARAM_DESCRIPTION, Default: TRUE
#' @param right_inclusive PARAM_DESCRIPTION, Default: FALSE
#'
#' @return OUTPUT_DESCRIPTION
#'
#' @details DETAILS
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#'
#' @seealso
#'  \code{\link[dplyr]{mutate}}
#'  \code{\link[rlang]{sym}}
#'
#' @rdname tile_metrics
#'
#' @importFrom dplyr mutate
#' @importFrom rlang sym
tile_metrics = function(metrics_df,
                        genome_tile_df, chr_seqlength,
                        cut_col = "POS",
                        left_inclusive = TRUE,
                        right_inclusive = FALSE){

  metrics_df %>%
    mutate(tile =
             cut(!!rlang::sym(cut_col),
                 c(genome_tile_df$start, chr_seqlength + 1),
                 include.lowest = left_inclusive,
                 right = right_inclusive))

}
