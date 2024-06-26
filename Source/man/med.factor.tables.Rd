% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ShedsHT.R
\name{med.factor.tables}
\alias{med.factor.tables}
\title{med.factor.tables}
\usage{
med.factor.tables(ef, media.sur)
}
\arguments{
\item{ef}{A data set created internally using the \code{\link{run}} function and the \code{\link{read.exp.factors}} function
to import the user specified Exp_factors input file. The data set contains the distributional parameters for the
exposure factors. All of these variables may have age or gender-dependent distributions, although in the absence of data, many
are assigned a single distribution from which all persons are sampled.}

\item{media.sur}{A list of surface media. This data set is created internally by sub-setting the Media input
file (read in with the \code{\link{read.media.file}} function within the \code{\link{run}} function) to extract only surface
media.}
}
\value{
exp.med Media specific exposure factors presented as a data table. The present version of the model consists of only
3 such exposure factors, but the code will accept more. Depending on the user input, these generated exposure factors may be
gender and/or season specific in addition to media specific. The output consists of 24 rows per variable (2 genders x 4
seasons x 3 media), when all ages share the same distribution. If a variable has N age categories (each with its own
distribution) then there are (24 x N) rows for that variable.
In addition, if the function runs successfully, the following message is printed: "Media-specific Factor Tables completed"
}
\description{
Constructs tables of media specific exposure factors for each relevant combination of age, gender, and season, and each
of three media specific variables in the \code{ef} argument. These are avail.f, dermal.tc, and om.ratio.
}
\details{
The three media specific variables in the \code{ef} argument are as follows:
\code{avail.f} Fraction of chemical available for transfer from surfaces via touching.
\code{dermal.tc} Dermal transfer coefficient, cm^2/hr.
\code{om.ratio} Ratio of object-to-mouth exposure to indirect dermal exposure.
}
\note{
The first input argument to \code{med.factor.tables} (ef) is created by reading in the Exp_factors
input file (specified on the Run file) with the \code{\link{read.exp.factors}} function within the
\code{\link{run}} function. The \code{media.sur} input is created by sub-setting the Media input file (read in
with the \code{\link{read.media.file}} function within the \code{\link{run}} function).
}
\seealso{
\code{\link{run}}, \code{\link{read.exp.factors}}, \code{\link{read.media.file}}
}
\author{
Kristin Isaacs, Graham Glen
}
