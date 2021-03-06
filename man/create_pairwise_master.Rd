% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scMeth_functions.R
\name{create_pairwise_master}
\alias{create_pairwise_master}
\title{Create pairwise comparisons}
\usage{
create_pairwise_master(cpg, digital = TRUE, cores_to_use = 2,
  calcdiff = TRUE)
}
\arguments{
\item{cpg}{A list of named data frames containing CpG calls. See details for required format of dataframes. Required.}

\item{digital}{Whether or not to discard non-binary CpG calls. Useful in single-cells as it's very unlikely that a single-cell contains a heterozygous methylation call. Defaults to TRUE.}

\item{calcdiff}{Whether or not to directly calculate the average difference (if TRUE), or to return a list of dataframes containing pairwise common CpGs (if FALSE). Defaults to TRUE. Usually you don't want false unless you wish to do something else with all the pairwise data}

\item{ncores}{Number of cores to parallelize over. Defaults to 1}
}
\value{
A list of dataframes if \code{calcdiff} is FALSE. Otherwise, a dataframe containing the pairwise dissimiarlties if \code{calcdiff} is TRUE.
}
\description{
Create pairwise comparisons between single-cells
}
\details{
Each dataframe containing CpG calls must have the following four columns:
1. Chromsome column, named "chr"
2. Start/Position column, named "start"
3. Percentage or fractional methylation column, named "meth" (between 0-100 or 0-1)
}

