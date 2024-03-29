% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{.getClusterDatatable}
\alias{.getClusterDatatable}
\title{Map each geneset to the cluster it belongs}
\usage{
.getClusterDatatable(cluster, gs_names, gs_description)
}
\arguments{
\item{cluster}{A \code{list} of clusters}

\item{gs_names}{A vector of geneset names}

\item{gs_description}{A vector of descriptions for each geneset}
}
\value{
A \code{data.frame} mapping each geneset to the cluster(s) it belongs to
}
\description{
Map each geneset to the cluster it belongs and return the information as
a \code{data.frame}
}
