% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parsers.R
\name{parse_magma_geneset}
\alias{parse_magma_geneset}
\title{Read in the results of magma geneset into tidy format}
\usage{
parse_magma_geneset(path)
}
\arguments{
\item{path}{filepath for gsa.out file from magma-geneset analysis}
}
\value{
a tibble
}
\description{
Read in the results of magma geneset into tidy format
}
\examples{
\dontrun{
parse_magma_geneset("/magma_geneset/amygdala_excitatory.gsa.out")
}
}
