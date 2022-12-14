% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bin_variants.R
\name{bin_variants}
\alias{bin_variants}
\title{Bin Variants}
\usage{
bin_variants(
  variant_metrics_df,
  tiled_genome_df,
  chr_seqlengths,
  chr_colname,
  qvalue_lower_thres = 0.1,
  qvalue_upper_thres = 0.5
)
}
\arguments{
\item{variant_metrics_df}{PARAM_DESCRIPTION}

\item{tiled_genome_df}{PARAM_DESCRIPTION}

\item{chr_seqlengths}{PARAM_DESCRIPTION}

\item{chr_colname}{PARAM_DESCRIPTION}

\item{qvalue_lower_thres}{PARAM_DESCRIPTION, Default: 0.1}

\item{qvalue_upper_thres}{PARAM_DESCRIPTION, Default: 0.5}
}
\value{
OUTPUT_DESCRIPTION
}
\description{
Calculate statistics over bins (genome tiles)
}
\details{
DETAILS
}
\examples{

\dontrun{
if(interactive()){
 #EXAMPLE1
 }
}
}
\seealso{
\code{\link[dplyr]{group_by}},
\code{\link[dplyr]{group_split}},
\code{\link[dplyr]{pull}},
\code{\link[dplyr]{summarise}},
\code{\link[dplyr]{context}},
\code{\link[dplyr]{mutate}},
\code{\link[dplyr]{select}},
\code{\link[dplyr]{reexports}}
\code{\link[rlang]{sym}}
\code{\link[purrr]{map}}
\code{\link[stringr]{str_remove}}
\code{\link[tidyr]{separate}}
}
