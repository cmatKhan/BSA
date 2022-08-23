# TODO THIS WAS USING BINNER, NOT BINNER2


#' @title plot bin dist
#' @description This function makes the distribution plots for the bins, using the output table of the binner()palette() function as input.
#'
#' @param inputBinTable PARAM_DESCRIPTION
#' @param windowsSizes PARAM_DESCRIPTION, Default: c(5000, 10000, 20000, 50000)
#' @param fileName PARAM_DESCRIPTION
#'
#' @return OUTPUT_DESCRIPTION
#'
#' @details DETAILS
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_density scale_x_continuous scale_color_manual theme element_text xlab ggsave
#' @importFrom grDevices palette
binDistPlotter <- function(inputBinTable, windowsSizes = c(5000, 10000, 20000, 50000), fileName) {

  binnerDataframe <- data.frame(
    CHR = as.character(), bin_floor = as.numeric(),
    bin_ceiling = as.numeric(), bin_middle = as.numeric(),
    lungAverage = as.numeric(), YPDaverage = as.numeric(),
    mut_perBin = as.numeric(), bin.size = as.numeric()
  )
  for (z in 1:length(windowsSizes)) {
    binnerDataframe <- rbind(binnerDataframe, binner2(inputBinTable, bin.size = windowsSizes[z]))
  }
  zeta <- melt(binnerDataframe, id.vars = colnames(binnerDataframe[, c(1:4, 7, 8)]))
  zeta_lung <- subset(zeta, variable == "lungMean")
  zeta_ypd <- subset(zeta, variable == "YPDmean")

  ggplot2::ggplot(data = zeta_lung, ggplot2::aes(x = (zeta_lung$value * 100),
                                                 color = factor(zeta_lung$Bin_size))) +
    ggplot2::geom_density() +
    # facet_wrap(~factor(zeta_lung$Chr,levels=unique(zeta_lung$Chr)), ncol = 3)+
    ggplot2::scale_x_continuous(limits = c(-100, 100)) +
    ggplot2::scale_color_manual(name = "Bin size",
                                labels = as.character(windowsSizes),
                                values = grDevices::palette()[1:length(windowsSizes)]) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 12),
      axis.title.y = ggplot2::element_text(size = 12),
      axis.text.x = ggplot2::element_text(size = 8),
      axis.text.y = ggplot2::element_text(size = 6),
      strip.text.x = ggplot2::element_text(size = 10)
    ) +
    ggplot2::xlab("Enrichment (Percentage points)")
  ggplot2::ggsave(paste0(fileName, ".binLungDistribution.png"))

  ggplot2::ggplot(data = zeta_ypd,
                  ggplot2::aes(x = (zeta_ypd$value * 100),
                               color = factor(zeta_ypd$Bin_size))) +
    ggplot2::geom_density() +
    # facet_wrap(~factor(zeta_lung$Chr,levels=unique(zeta_lung$Chr)), ncol = 3)+
    ggplot2::scale_x_continuous(limits = c(-100, 100)) +
    ggplot2::scale_color_manual(name = "Bin size",
                                labels = as.character(windowsSizes),
                                values = grDevices::palette()[1:length(windowsSizes)]) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 12),
      axis.title.y = ggplot2::element_text(size = 12),
      axis.text.x = ggplot2::element_text(size = 8),
      axis.text.y = ggplot2::element_text(size = 6),
      strip.text.x = ggplot2::element_text(size = 10)
    ) +
    ggplot2::xlab("Enrichment (Percentage points)")
  ggplot2::ggsave(paste0(fileName, ".binYpdDistribution.png"))
}
