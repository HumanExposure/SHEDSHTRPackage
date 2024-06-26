% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utility.R
\name{summary.stats}
\alias{summary.stats}
\title{summary.stats}
\usage{
summary.stats(x)
}
\arguments{
\item{x.}{Data set passed from summarize.chemical for a cohort}
}
\value{
y A data frame object with 23 rows and 11 columns, with each column being an exposure or dose variable,
and each row containing a statistic for that variable. For each expsoure varianle the total expsoure, quantiles, mean, and
SD.
}
\description{
Summary.stats constructs the table of exposure and dose statistics for cohort entered into \code{\link{summarize.chemical}}.
}
\details{
Summary.stats is called by \code{\link{summarize.chemical}}.  The input data set "y"
is one population cohort from the exposure data set passed into \code{\link{summarize.chemical}}.
}
\seealso{
\code{\link{summarize.chemical}}
}
\author{
Kristin Isaacs, Graham Glen
}
