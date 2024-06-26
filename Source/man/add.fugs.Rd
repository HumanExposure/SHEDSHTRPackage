% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Fugacity.R
\name{add.fugs}
\alias{add.fugs}
\title{add.fugs}
\usage{
add.fugs(n.per, x, pdmf)
}
\arguments{
\item{n.per}{The total number of simulated persons in this model run specified in the Run file}

\item{x}{the output of the \code{\link{read.fug.inputs}} function.}

\item{pdmf}{he output of the \code{\link{add.factors}} function. A data set containing physiological and demographic parameters for each theoretical person, the duration of exposure to each potential exposure medium for each person,
the media specific exposure factors, and the  number of baths taken and hand wash events occurring per day for each person.}
}
\value{
pdmff Output contains values sampled from the distributions of each relevant variable in the Fugacity
 input file for each theoretical person.
}
\description{
Evaluates the variables in the fugacity input file and creates values for those variables corresponding to each simulated person.
}
\details{
This function evaluates the variables in the fugacity input file, creating one value for each simulated
person and adding each variable as a new column in the pdmf data set (renamed as pdmff). All variables are left i
n their original units except for those with units ug/cm2, which are converted to ug/m2. In the fugacity calculations,
all masses are in ug and all lengths are in m.
}
\seealso{
\code{\link{add.factors}}, \code{\link{add.media}}, \code{\link{run}}
}
\author{
Kristin Isaacs, Graham Glen
}
\keyword{SHEDS}
