% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{.checkPPI}
\alias{.checkPPI}
\title{Check PPI format}
\usage{
.checkPPI(ppi)
}
\arguments{
\item{ppi}{a \code{data.frame}, Protein-protein interaction (PPI) network data frame.
The object is expected to have three columns, \code{Gene1} and \code{Gene2}
which specify the gene names of the interacting proteins in no
particular order (symmetric interaction) and a column
\code{combined_score} which is a numerical value of the strength of
the interaction.}
}
\value{
A validated and formatted PPI data frame.
}
\description{
Check if the Protein-Protein-interaction (PPI) has the expected format for
this app
}
