% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ShedsHT.R
\name{create.scen.factors}
\alias{create.scen.factors}
\title{create.scen.factors}
\usage{
create.scen.factors(f)
}
\arguments{
\item{f}{An internally generated data set from the Source_vars input file (as specified in the Run
file) containing data on the distribution of each source variable in SHEDS.HT.}
}
\value{
dt A data set with values and probabilities associated with each variable distribution presented in \code{f}.
All variables in \code{f} are also retained.
}
\description{
This function takes the information on distributions from the \code{all.scenarios} data set (which comes from
the Source_variables input file) and converts it into the parameter set needed by SHEDS.HT.
}
\details{
The steps involved in this function are 1) converting \code{prevalence} from a percentage to a binomial form, 2)
converting \code{CV} to standard deviation for normals, and 3) converting \code{Mean} and \code{CV} to par1 and par2 for
lognormals.
Note that \code{prevalence} in SHEDS becomes a binomial distribution, which returns a value of either 0 or 1 when evaluated.
Each simulated person either "does" or "does not" partake in this scenario. Similar logic applies to the \code{frequency}
variable, except that the returned values may be larger than one (that is, 2 or more) for very frequent scenarios.  All the
exposure equations contain the \code{prevalence} variable. If \code{prevalence} is set to one for that person, the exposure
is as expected, but if \code{prevalence}=0 then no exposure occurs.
}
\author{
Kristin Isaacs, Graham Glen
}
