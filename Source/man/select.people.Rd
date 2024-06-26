% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ShedsHT.R
\name{select.people}
\alias{select.people}
\title{select.people}
\usage{
select.people(n, pop, py, act.p, diet.p, act.d, diet.d,specs)
}
\arguments{
\item{n}{Number of persons.}

\item{pop}{The population input; output of \code{\link{read.pop.file}} function. Contains counts by gender and each year of age
from the 2000 U.S. census.  When a large age range is modeled, this ensures that SHEDS chooses age and gender with the correct
overall probability.
#'}

\item{py}{Regression parameters on the three physiological variables of interest (weight, height, and body mass index) for
various age and gender groups. Output of the \code{\link{read.phys.file}} function.}

\item{act.p}{Activity diary pools; output of the \code{\link{act.diary.pools}} function. Each element in this input is a list
of acceptable activity diary numbers for each year of age, for each gender, weekend, and season combination.}

\item{diet.p}{Dietary diary pools; output of the \code{\link{diet.diary.pools}} function. Each element in this input is a list
of acceptable activity diary numbers for each year of age, for each gender, weekend, and season combination.}

\item{act.d}{Activity diaries, which indicate the amount of time and level of metabolic activity in various 'micros'.
Each line of data represents one person-day (24 hours). Output of the \code{\link{read.act.diaries}} function.}

\item{diet.d}{Daily diaries of dietary consumption by food group; output of the \code{\link{read.diet.diaries}} function.
Each line represents one person-day, with demographic variables followed by amounts (in grams/day) for a list of food types
indicated by a short abbreviation on the header line.}

\item{specs}{Output of the \code{\link{read.run.file}} function, which can be modified by the \code{\link{update.specs}}
function before input into \code{\link{select.people}}.}
}
\value{
pd A dataframe of "person-demographics": assigned demographic and physiological parameters for each theoretical
person modeled in SHEDS.HT
}
\description{
Assigns demographic and physiological variables to each theoretical person to be modeled.
}
\details{
This is the first real step in the modeling process. It first fills an array \code{q} with uniform random numbers,
with ten columns because there are 10 random variables defined by this function. There number of rows correspond to the number
of persons, capped at the \code{set.size} specified in the Run input file (typically 5000).  Gender is selected
from a discrete (binomial) distribution where the counts of males and females in the study age range determines the gender
probabilities. Age is tabulated next, separately for each gender. The counts by year of age are chosen for the appropriate
gender and used as selection weights. Season is assigned randomly (equal weights) using those specified in the Run
input file. \code{Weekend} is set to one or zero, with a chance of 2/7 for the former.
The next block of code assigns physiological variables. Weight is lognormal in SHEDS, so a normal is sampled first and then
\code{exp()} is applied.  This means that the weight parameters refer to the properties of log(weight), which were fit by
linear regression. The basal metabolic rate (bmr), is calculated by regression. A minimum bmr is set to prevent extreme cases
from becoming zero or negative.  The alveolar breathing ventilation rate corresponding to bmr is also calculated. The SHEDS
logic sets activities in each micro to be a multiple of these rates, with outdoor rates higher than indoor, and indoor rates
higher than sleep rates. This calculated activities affect the inhaled dose. The skin surface area is calculated using
regressions based on height and weight for 3 age ranges.
The next step is to assign diaries.  Here, a FOR loop over n persons (n rows) is used to assign appropriate diet and activity
diary pools to each person. An empirical distribution is created,  consisting of the list of diary numbers for each pool.
The final step is to retrieve the actual data from the chosen activity and diet diaries, and the result becomes \code{pd}.
}
\seealso{
\code{\link{run}}, \code{\link{read.pop.file}}, \code{\link{read.phys.file}}, \code{\link{act.diary.pools}}, \code{\link{diet.diary.pools}}, \code{\link{read.act.diaries}}, \code{\link{read.diet.diaries}}, \code{\link{update.specs}}
}
\author{
Kristin Isaacs, Graham Glen
}
