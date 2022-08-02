#' @title a simulate function
#' @description FUNCTION_DESCRIPTION
#'
#' @param SNPset PARAM_DESCRIPTION
#' @param popStruc PARAM_DESCRIPTION, Default: 'F2'
#' @param bulkSize PARAM_DESCRIPTION
#' @param depth PARAM_DESCRIPTION, Default: 1:100
#' @param replications PARAM_DESCRIPTION, Default: 10000
#' @param filter PARAM_DESCRIPTION, Default: 0.3
#' @param intervals PARAM_DESCRIPTION, Default: c(0.05, 0.025)
#'
#' @return OUTPUT_DESCRIPTION
#'
#' @details DETAILS
#'
#' @export
#'
#' @importFrom timeSeries quantile as.data.frame t
simulateConfInt <-
  function(SNPset, popStruc = "F2",
           bulkSize,
           depth = 1:100,
           replications = 10000,
           filter = 0.3,
           intervals = c(0.05, 0.025)) {
    if (popStruc == "F2") {
      message(
        "Assuming bulks selected from F2 population, with ",
        as.character(bulkSize),
        " individuals per bulk."
      )
    } else {
      message(
        "Assuming bulks selected from RIL population, with ",
        as.character(bulkSize),
        " individuals per bulk."
      )
    }

    # makes a vector of possible alt allele frequencies once. this is then sampled for each replicate

    # This part commented below was replaced in the call below.

    # tmp_freq <-
    # replicate(n = replications * 10, simulateAlleleFreq(n = bulkSize, pop = popStruc))

    message(
      paste0(
        "Simulating ",
        as.character(replications),
        " SNPs with reads at each depth: ",
        as.character(min(depth)),
        "-",
        as.character(max(depth))
      )
    )
    message(paste0(
      "Keeping SNPs with >= ",
      as.character(filter),
      " SNP-index in both simulated bulks"
    ))

    # tmp allele freqs are sampled to produce 'replicate' numbers of probablities. these
    # are then used as altFreq probs to simulate SNP index values, per bulk.
    CI <- sapply(
      X = depth,
      FUN = function(x) {
        timeSeries::quantile(
          x = BSA::simulateSNPindex(
            depth = x,
            altFreq1 = SNPset$SNPindex.LOW, # I adapted this. Instead of sampling the tmp_freq, we are using the SNP index from the inoculum to calculate this.

            altFreq2 = SNPset$SNPindex.LOW,
            replicates = replications,
            filter = filter
          ),
          probs = intervals,
          names = TRUE
        )
      }
    )

    CI <- timeSeries::as.data.frame(CI)

    if (length(CI) > 1) {
      CI <- data.frame(timeSeries::t(CI))
    }

    names(CI) <- paste0("CI_", 100 - (intervals * 200))
    CI <- cbind(depth, CI)

    # to long format for easy plotting
    # tidyr::gather(data = CI,
    #     key = interval,
    #     convert = TRUE,
    #     value = SNPindex,-depth) %>%
    #     dplyr::mutate(Confidence = factor(ifelse(
    #         interval > 0.5,
    #         paste0(round((1 - interval) * 200, digits = 1), "%"),
    #         paste0((interval * 200), "%")
    # )))
    CI
  }
