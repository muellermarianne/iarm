% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Person-Fit.R
\name{test_prop}
\alias{test_prop}
\title{Properties of the Test}
\usage{
test_prop(object)
}
\arguments{
\item{object}{An object of class "Rm", a fitted Rasch model or partial
credit model using  the functions RM or PCM in package eRm, or an object of class "pcmodel",
 a fitted partial credit model using the function pcmodel in package psychotools.}
}
\value{
a list containing:
\item{Separation reliability}{the person separation reliability as calculated in package eRm for objects of class "Rm".}
\item{Test difficulty}{person value with an expected score equal to half of the maximum score.}
\item{Test target}{person value where test information is maximized.}
\item{Test information}{maximal value of the test information}
}
\description{
Information summarizing measurement quality of the test and test targeting.
}
\examples{
rm.mod <- RM(amts[,4:13])
test_prop(rm.mod)
}
\references{
Christensen, K. B. , Kreiner, S. & Mesbah, M. (Eds.)
\emph{Rasch Models in Health}. Iste and Wiley (2013), pp. 63 - 70.
}
\author{
Marianne Mueller
}
