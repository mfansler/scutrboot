
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scUTRboot

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/mfansler/scutrboot)](https://github.com/mfansler/scutrboot/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/mfansler/scutrboot)](https://github.com/mfansler/scutrboot/pulls)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8057843.svg)](https://doi.org/10.5281/zenodo.8057843)

<!-- badges: end -->

## Overview

The `scutrboot` package provides a collection of statistical tools for
the analysis of 3’ UTR isoform usage from scRNA-seq data. In particular,
it was developed to accompany the scUTRquant pipeline, working from 3’
UTR isoform-level count data to characterize isoform usage in groups of
cells (e.g., clusters, cell types) and test differential usage across
group. The “*boot*”, both refers to the package’s role as the part of
the scUTRx suite where the tools are kept and to the bootstrap-based
procedures that are core to the statistical framework.

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `scutrboot` from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# not available yet
# BiocManager::install("scutrboot")
```

And the development version from
[GitHub](https://github.com/mfansler/scutrboot) with:

``` r
BiocManager::install("mfansler/scutrboot")
```

## Usage

Application of the `scutrboot` package can be found in [the
`scUTRquant-figures`
repository](https://github.com/Mayrlab/scUTRquant-figures), particularly
in the notebooks used to generate Figures
[4](https://github.com/Mayrlab/scUTRquant-figures/tree/main/figures/figure4)
and
[5](https://github.com/Mayrlab/scUTRquant-figures/tree/main/figures/figure5).

## Citation

Below is the citation output from using `citation('scutrboot')` in R.
Please run this yourself to check for any updates on how to cite
**scutrboot**.

``` r
print(citation('scutrboot'), bibtex = TRUE)
#> 
#> To cite package 'scutrboot' in publications use:
#> 
#>   Mervin Fansler (2021). scutrboot: Single-Cell UTR Bootstrap Tools. R
#>   package version 0.2.2.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {scutrboot: Single-Cell UTR Bootstrap Tools},
#>     author = {Mervin Fansler},
#>     year = {2021},
#>     note = {R package version 0.2.2},
#>   }
```

Please note that the `scutrboot` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `scutrboot` project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

-   Continuous code testing is possible thanks to [GitHub
    actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
    through *[usethis](https://CRAN.R-project.org/package=usethis)*,
    *[remotes](https://CRAN.R-project.org/package=remotes)*, and
    *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)*
    customized to use [Bioconductor’s docker
    containers](https://www.bioconductor.org/help/docker/) and
    *[BiocCheck](https://bioconductor.org/packages/3.14/BiocCheck)*.
-   Code coverage assessment is possible thanks to
    [codecov](https://codecov.io/gh) and
    *[covr](https://CRAN.R-project.org/package=covr)*.
-   The [documentation website](http://mfansler.github.io/scutrboot) is
    automatically updated thanks to
    *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
-   The code is styled automatically thanks to
    *[styler](https://CRAN.R-project.org/package=styler)*.
-   The documentation is formatted thanks to
    *[devtools](https://CRAN.R-project.org/package=devtools)* and
    *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.14/biocthis)*.
