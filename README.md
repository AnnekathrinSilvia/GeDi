
<img src="man/figures/GeDi.png" align="right" alt="" width="120" />

<!-- README.md is generated from README.Rmd. Please edit that file -->

# GeDi

<!-- badges: start -->

[![R build
status](https://github.com/AnnekathrinSilvia/GeDi/workflows/R-CMD-check/badge.svg)](https://github.com/AnnekathrinSilvia/GeDi/actions)

[![](https://img.shields.io/github/last-commit/AnnekathrinSilvia/GeDi.svg)](https://github.com/AnnekathrinSilvia/GeDi/commits/devel)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![CodeCov.io codecov
status](https://codecov.io/github/AnnekathrinSilvia/GeDi/coverage.svg?branch=master)](https://codecov.io/github/AnnekathrinSilvia/GeDi)
<!-- badges: end -->

The goal of GeDi is to provide a user-friendly platform for exploring
and analyzing the results of functional annotation and enrichment
analyses. By integrating various enrichment methods and offering
interactive visualizations, GeDi aims to streamline the process of
interpreting biological data and identifying meaningful patterns within
it. Through features like distance score calculation, clustering, and
bookmarking, GeDi empowers users to gain deeper insights into their data
and make informed decisions in their research endeavors.

## Installation

You can install the development version of GeDi from
[GitHub](https://github.com/AnnekathrinSilvia/GeDi) with:

``` r
install.packages("devtools")
devtools::install_github("AnnekathrinSilvia/GeDi")
```

## Example

If you want to give GeDi a testrun on a demo dataset, you can simply
launch the Shiny application and load the demo data in the **Data
input** panel. To start the Shiny Application, you can use the following
code:

``` r
library("GeDi")
GeDi()
```

## Usage Overview

You can find the rendered version of the documentation of `GeDi` at the
project website <https://AnnekathrinSilvia.github.io/GeDi>, created with
`pkgdown`.

## Development

If you encounter a bug, have usage questions, or want to share ideas and
functionality to make this package better, feel free to file an
[issue](https://github.com/AnnekathrinSilvia/GeDi/issues).

## Code of Conduct

Please note that the GeDi project is released with a [Contributor Code
of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## License

MIT © Annekathrin Silvia Nedwed
