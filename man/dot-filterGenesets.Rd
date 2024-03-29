% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{.filterGenesets}
\alias{.filterGenesets}
\title{Filter Genesets from the input data}
\usage{
.filterGenesets(remove, df_genesets)
}
\arguments{
\item{remove}{a \code{list}, A list of geneset names to be removed}

\item{df_genesets}{a \code{data.frame}, A \code{data.frame} with at least two columns.
One should be called \code{Geneset}, containing the
names/identifiers of the genesets in the data. The second
column should be called \code{Genes} and contains one string of
the genes contained in each geneset.}
}
\value{
A \code{data.frame} containing information about filtered genesets
}
\description{
Filter a preselected list of genesets from a \code{data.frame} of genesets
}
