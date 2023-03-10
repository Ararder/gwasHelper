% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cleansumstats.R
\name{migrate}
\alias{migrate}
\title{Title}
\usage{
migrate(infile, meta, name)
}
\arguments{
\item{infile}{raw gwas sumstat to clean}

\item{meta}{filepath to metafile}

\item{name}{gwas sumstat name}
}
\value{
a filepath to the slurm job
}
\description{
Title
}
\examples{
\dontrun{
migrate("my_raw_gwas.tsv.gz",name =  "mdd2024",col_map = col_map)
}
}
