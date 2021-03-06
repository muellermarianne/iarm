% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Item-Fit.R
\name{boot_fit}
\alias{boot_fit}
\title{Computes Bootstrapping P Values for Outfit and Infit Statistics}
\usage{
boot_fit(
  object,
  B,
  p.adj = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "none")
)
}
\arguments{
\item{object}{an object of class "Rm" (output of RM or PCM) or class "pcmodel"}

\item{B}{Number of replications.}

\item{p.adj}{Correction method for multiple testing. The methods are "BH","holm", "hochberg", "hommel", "bonferroni", "BY", "none". See \code{\link{p.adjust}}.}
}
\value{
object of class bootfit with outfit and infit statistics and corresponding p values.
}
\description{
Computes Bootstrapping P Values for Outfit and Infit Statistics
}
