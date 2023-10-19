
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BSA

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/cmatKhan/BSA)](https://github.com/cmatKhan/BSA/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/cmatKhan/BSA)](https://github.com/cmatKhan/BSA/pulls)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/cmatKhan/BSA/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/cmatKhan/BSA/actions/workflows/check-bioc.yml)
[![Codecov test
coverage](https://codecov.io/gh/cmatKhan/BSA/branch/master/graph/badge.svg)](https://app.codecov.io/gh/cmatKhan/BSA?branch=master)
[![Software
DOI](https://zenodo.org/badge/520266986.svg)](https://zenodo.org/badge/latestdoi/520266986)
<!-- badges: end -->

The goal of `BSA` is to reimplement Daniel’s BSA2 code in a package
structure such that it may be installed using R `install`. This may also
form the seed of either a software or
[workflow](https://contributions.bioconductor.org/non-software.html)
package on bioconductor.

## Installation

``` r

install.packages("remotes")

remotes::install_github("cmatKhan/BSA")
```

## Usage

The first thing that you will need to do iscreate a samplesheet and
process the data through the [BSA processing
pipeline](https://github.com/BrentLab/callvariants).

Instructions for creating this samplesheet are in the Article (above)
called `CreatePipelineSampleSheet`.

Next, you will be using the functions in this package to analyze the
data. An example of the process are in `BSA3` and `BSA6`.

## Citation

Below is the citation output from using `citation('BSA')` in R. Please
run this yourself to check for any updates on how to cite **BSA**.

``` r
print(citation('BSA'), bibtex = TRUE)
#> To cite package 'BSA' in publications use:
#> 
#>   Agustihno, Daniel (2023). "Daniel's paper placeholder." _bioRxiv_. ,
#>   <https://www.biorxiv.org/content/10.1101/TODO>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Daniel's paper placeholder},
#>     author = {{Agustihno} and {Daniel}},
#>     year = {2023},
#>     journal = {bioRxiv},
#>     doi = {10.1101/TODO},
#>     url = {https://www.biorxiv.org/content/10.1101/TODO},
#>   }
```

Please note that the `BSA` was only made possible thanks to many other R
and bioinformatics software authors, which are cited either in the
vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `BSA` project is released with a [Contributor Code
of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

- Continuous code testing is possible thanks to [GitHub
  actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
  through *[usethis](https://CRAN.R-project.org/package=usethis)*,
  *[remotes](https://CRAN.R-project.org/package=remotes)*, and
  *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized
  to use [Bioconductor’s docker
  containers](https://www.bioconductor.org/help/docker/) and
  *[BiocCheck](https://bioconductor.org/packages/3.17/BiocCheck)*.
- Code coverage assessment is possible thanks to
  [codecov](https://codecov.io/gh) and
  *[covr](https://CRAN.R-project.org/package=covr)*.
- The [documentation website](http://cmatKhan.github.io/BSA) is
  automatically updated thanks to
  *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
- The code is styled automatically thanks to
  *[styler](https://CRAN.R-project.org/package=styler)*.
- The documentation is formatted thanks to
  *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.17/biocthis)*.
