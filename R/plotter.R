#TODO complete documentation

#' @title plot some stuff
#' @description This function will plot the enrichment of each SNP.
#'   It uses the output of the percenter() as input.
#'   Obs, I have to set to run the percenter() as well as define the
#'   fileNames outside the function, in the for loop to make it work.
#'   here what I mean:  l=percenter(sampleColumn=sampleColumn)
#'
#' @param input PARAM_DESCRIPTION
#' @param fileName PARAM_DESCRIPTION
#'
#' @return OUTPUT_DESCRIPTION
#'
#' @details DETAILS
#'
#' @export
#' @importFrom ggplot2 ggplot geom_point aes scale_color_manual facet_wrap geom_hline scale_y_continuous guides guide_legend xlab ylab theme element_text ggsave
plotter <- function(input, fileName) {
  l <- input

  ggplot2::ggplot(data = l) +
    ggplot2::geom_point(ggplot2::aes(x = POS, y = (Lung * 100),
                                     colour = "blue"), size = 0.01) +
    ggplot2::geom_point(ggplot2::aes(x = POS, y = (YPD * 100),
                                     colour = "orange"), size = 0.01) +
    ggplot2::scale_color_manual(name = "Condition",
                                labels = c("Lung", "YPD"),
                                values = c("blue", "orange")) +
    ggplot2::facet_wrap(~ factor(l$CHR, levels = unique(l$CHR)),
                        scales = "free", ncol = 3) +
    ggplot2::geom_hline(data = l, ggplot2::aes(yintercept = 0)) +
    ggplot2::scale_y_continuous(limits = c(-100, 100)) +
    ggplot2::guides(colour =
                      ggplot2::guide_legend(override.aes = list(size = 6))) +
    ggplot2::xlab("Position") +
    ggplot2::ylab("Enrichment (Percentage points)") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 18),
      axis.title.y = ggplot2::element_text(size = 18),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.text.y = ggplot2::element_text(size = 12),
      strip.text.x = ggplot2::element_text(size = 14)
    )
  ggplot2::ggsave(paste0(fileName, ".png"),
                  width = 18, units = "in", height = 14)

  ggplot2::ggplot(data = l) +
    ggplot2::geom_point(ggplot2::aes(x = POS,
                                     y = (Lung * 100),
                                     colour = "blue"), size = 0.001) +
    ggplot2::scale_color_manual(name = "Condition",
                                labels = c("Lung"), values = c("blue")) +
    ggplot2::facet_wrap(~ factor(l$CHR,
                                 levels = unique(l$CHR)),
                        scales = "free", ncol = 3) +
    ggplot2::geom_hline(data = l, ggplot2::aes(yintercept = 0)) +
    ggplot2::scale_y_continuous(limits = c(-100, 100)) +
    ggplot2::xlab("Position") +
    ggplot2::ylab("Enrichment (Percentage points)") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 18),
      axis.title.y = ggplot2::element_text(size = 18),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.text.y = ggplot2::element_text(size = 12),
      strip.text.x = ggplot2::element_text(size = 14)
    )
  ggplot2::ggsave(paste0(fileName, ".lung_only.png"),
                  width = 18, units = "in", height = 14)

  ggplot2::ggplot(data = l) +
    ggplot2::geom_point(ggplot2::aes(x = POS, y = (YPD * 100),
                                     colour = "orange"), size = 0.001) +
    ggplot2::scale_color_manual(name = "Condition",
                                labels = c("YPD"), values = c("Orange")) +
    ggplot2::facet_wrap(~ factor(l$CHR, levels = unique(l$CHR)),
                        scales = "free", ncol = 3) +
    ggplot2::geom_hline(data = l, ggplot2::aes(yintercept = 0)) +
    ggplot2::scale_y_continuous(limits = c(-100, 100)) +
    ggplot2::xlab("Position") +
    ggplot2::ylab("Enrichment (Percentage points)") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 18),
      axis.title.y = ggplot2::element_text(size = 18),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.text.y = ggplot2::element_text(size = 12),
      strip.text.x = ggplot2::element_text(size = 14)
    )
  ggplot2::ggsave(paste0(fileName, ".YPD_only.png"),
                  width = 18, units = "in", height = 14)
}
