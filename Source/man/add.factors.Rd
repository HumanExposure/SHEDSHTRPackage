% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ShedsHT.R
\name{add.factors}
\alias{add.factors}
\title{add.factors
Adds specific exposure factors to the \code{pdm} data table, which is output from the \code{\link{add.media}} function. Here,
"specific" means taking into account the age, gender, season, and exposure media for each person.}
\usage{
add.factors(n, gen.f, med.f, exp.f, surf, pdm)
}
\arguments{
\item{n}{Number of persons}

\item{gen.f}{Non-media specific exposure factors as a data table. Output from the \code{\link{gen.factor.tables}} function.}

\item{med.f}{Media specific exposure factors presented as a data table. Output from the \code{\link{med.factor.tables}} function.}

\item{exp.f}{Distributional parameters for the exposure factors. Output of the exp.factors function.}

\item{surf}{A list of surface media. Modified output of the \code{\link{read.media.file}} function.}

\item{pdm}{A data table containing the \code{pd} data frame of physiological and demographic parameters for each theoretical
person, and the \code{dur} array, which specifies the duration of exposure to each potential exposure medium for each person
in  \code{pd}. Output of the \code{\link{add.media}} function.}
}
\value{
pdmf A data set containing the \code{pdm} data table as well as media specific exposure factors, the number of baths
taken, and the number of hand wash events occurring per day per person contained in \code{pdm}.
}
\description{
add.factors
Adds specific exposure factors to the \code{pdm} data table, which is output from the \code{\link{add.media}} function. Here,
"specific" means taking into account the age, gender, season, and exposure media for each person.
}
\details{
The process of adding specific exposure factors to \code{pdm} involves multiple steps. First, \code{w} is determined,
which is the number of general factors plus the product of the number of media-specific factors and the number of surface media.
Air media do not have media specific factors in this version of SHEDS. An array, \code{q}, of uniform random samples is
generated, with one row per person and \code{w} columns.  A zero matrix, \code{r}, of the same size is defined.
Once these matrices are defined, the media-specific factors are determined. Two nested loops over variable and surface type
generate the values, which are stored in \code{r}. Next, another FOR loop determines the general factors. The \code{p} data set
contains the age, gender, and season for each person. These two data sets are then merged. The evaluation of these factors is
handled by the \code{\link{eval.factors}} function.
One of the exposure factors is \code{handwash.freq}. This was also part of SHEDS-Multimedia, where it represented the mean
number of hours in the day with hand washing events. An important aspect of that model was that because each person was followed
longitudinally, the actual number of hand washes on each day varied from one day to the next. Because of this, the distribution
for \code{handwash.freq} did not need to be restricted to integer values, as (for example) a mean of 4.5 per day is acceptable
and achievable, while choosing integer numbers of hand washes each day.  One of the early goals with SHEDS.HT was to attempt
to reproduce selected results from SHEDS-Multimedia. Therefore, similar logic was built into the current model. The
\code{hand.washes} variable is sampled from a distribution centered on \code{handwash.freq}, and then rounded to the nearest
integer.
The \code{bath} variable is another difficult concept. In theory, baths and showers are recorded on the activity diaries.
In practice, the activity diaries were constructed from approximately 20 separate studies, some of which did not contain enough
detail to identify separate bath or shower events. The result is that about half of all diaries record such events, but the
true rate in the population is higher. The \code{bath.p} variable was created to address this. It represents the  probability
that a non bath/shower activity diary should actually have one. Therefore, if the diary has one, then SHEDS automatically has
one. Otherwise, a binomial sample using \code{bath.p} as the probability is drawn. A bath/shower occurs unless both of these
are zero.
The effectiveness of hand washes or bath/shower at removing chemical from the skin is determined in the
\code{\link{post.exposure}} function.
}
\seealso{
\code{\link{eval.factors}}, \code{\link{post.exposure}}
}
\author{
Kristin Isaacs, Graham Glen
}
