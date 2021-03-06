% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Gamma.R
\name{partgam}
\alias{partgam}
\title{Conditional and Partial Gamma Coefficients}
\usage{
partgam(x, y, z, conf.level = 0.95)
}
\arguments{
\item{x, y, z}{Three numeric vectors or factors.}

\item{conf.level}{Confidence level for the returned confidence interval.}
}
\value{
data frame with estimates, standard errors and confidence interval limits.
}
\description{
Calculates conditional and partial Gamma coefficients for x and y given z with confidence intervals.
}
\references{
Davis, J. A. A Partial coefficient for Goodman and Kruskal's Gamma.
 \emph{Journal of the American Statistical Association}, 62 (317), 1967, pp. 189-193.
}
\seealso{
\code{\link{partgam_DIF}}, \code{\link{partgam_LD}}
}
\author{
Marianne Mueller
}
